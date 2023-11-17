#! /usr/bin/env python
"""
Replacement for JWST pipeline step for sky matching.

:Authors: Varun Bajaj


"""
import asdf
import logging
import numpy as np


from jwst.stpipe import Step

from jwst.datamodels import ModelContainer
from jwst.pipeline import calwebb_image3
from jwst.resample import resample_utils
from jwst.skymatch.skystatistics import SkyStats
from multiprocessing import Pool, cpu_count

from astropy.table import Table, vstack
from scipy.stats import sigmaclip

__all__ = ['PixelSkyMatchStep']


class PixelSkyMatchStep(Step):

    class_alias = "skymatch"

    spec = """
        # General sky matching parameters:
        skymethod = option('match', 'global+match', default='match') # sky computation method, more methods coming soon
        match_down = boolean(default=True) # adjust sky to lowest measured value?
        subtract = boolean(default=False) # subtract computed sky from image data?
        weight = boolean(default=True) # Weight the tile offsets by overlap?
        grouping = option('visit', None, default='visit')

        # Sky statistics parameters:
        skystat = option('median', 'midpt', 'mean', 'mode', default='median') # sky statistics
        lsigma = float(min=0.0, default=2.0) # Lower clipping limit, in sigma
        usigma = float(min=0.0, default=2.0) # Upper clipping limit, in sigma
        nclip = integer(min=0, default=5) # number of sky clipping iterations
        binwidth = float(min=0.0, default=0.1) # Bin width for 'mode' and 'midpt' `skystat`, in sigma

        # Other parameters
        save_tiles = boolean(default=False) # Save the drizzled images for each tile?
    """

    reference_file_types = []

    def process(self, input):
        self.log.setLevel(logging.DEBUG)
        if self.skymethod != 'match':
            raise NotImplementedError('Other sky methods are not implemented yet, but are coming soon.')

        # I don't understand what this does, but it's probably important
        # img = ModelContainer(
        #     input,
        #     save_open=not self._is_asn,
        #     return_open=not self._is_asn
        # )

        img = ModelContainer(input)

        # Get a name for this product
        filt = img[0].meta.instrument.filter
        pupil = img[0].meta.instrument.pupil
        det = img[0].meta.instrument.detector

        if hasattr(img, 'asn_table_name') and img.asn_table_name:
            aname = img.asn_table_name
            pname = '_'.join(aname.split('_')[:-1])
        elif hasattr(img.meta, 'filename') and img.meta.filename:
            aname = img.meta.filename
            pname = '_'.join(aname.split('_')[:-1])
        else:
            if pupil:
                pname = f'{filt}_{pupil}_{det}'
            else:
                pname = f'{filt}_{det}'

        # This creates a WCS/pixel grid that fits all the input exposures
        # This is needed to drizzle all the images to the same grid
        gw = resample_utils.make_output_wcs(img)
        tree = {"wcs": gw}

        wcs_file = asdf.AsdfFile(tree)
        # This seems unncecssary, perhaps there's a better way to specify the WCS params
        # The file must be written to disk for the resample step to use it
        self.wcs_filename = f'{pname}_dummy_wcs.asdf'
        fo = open(self.wcs_filename, 'wb')
        wcs_file.write_to(fo)
        fo.close()

        if self.grouping == 'visit':
            mcs = self.group_exps_by_visit(img)
        elif self.grouping is None:
            mcs = [ModelContainer([im]) for im in img]

        # Drizzle each tile and store the results

        for i, mc in enumerate(mcs):
            mc.meta.filename = f'tile{i}_{pname}'
        tiles = [self.driz_tile(mc) for mc in mcs]

        # No dreams of multiprocessing with Pool unless you can pickle datamodels
        # if self.n_process == 1:
        #     tiles = [self.driz_tile(mc) for mc in mcs]
        # else:
        #     if self.n_process == 0:
        #         n_cpu = min(cpu_count(), len(mcs))
        #     else:
        #         n_cpu = self.n_process
        #
        #     with Pool(n_cpu) as p:
        #         p.map(self.driz_tile, mcs)

        self._skystat = SkyStats(
            skystat=self.skystat,
            nclip=self.nclip,
            lsig=self.lsigma,
            usig=self.usigma,
            binwidth=self.binwidth
        )

        # Create matrices and perform least squares regression to find tile backgrounds
        equations, pdiffs, overlaps = self._create_matrix(tiles)
        if self.weight:
            W = np.sqrt(np.diag(overlaps/np.mean(overlaps)))
            equations = np.dot(W, equations)
            pdiffs = np.dot(pdiffs, W)
        lres = np.linalg.lstsq(equations, pdiffs, rcond=1.0E-12)
        # print(lres)
        tile_bgs = lres[0]

        # Adjust background levels to start from zero
        if self.match_down:
            matched_bgs = tile_bgs - np.amin(tile_bgs)
        else:
            ## TODO: Test matching up
            matched_bgs = tile_bgs - np.amax(tile_bgs)

        # Combine background levels with exposure offsets
        tbls = []
        for tile, bg in zip(tiles, matched_bgs):
            tbls.append(self.combine_with_exp(tile, bg))

        # Stack tables and write the results to a file
        bkgtbl = vstack(tbls)
        bkgtbl.write(f'{pname}_skies.txt', format='ascii.commented_header', overwrite=True)

        bg_dict = {row['FILENAME']:row['MATCHEDBG'] for row in bkgtbl}
        print(bkgtbl)

        self._set_sky_levels(img, bg_dict)

        return img

    def driz_tile(self, mc):
        """
        Drizzle exposures for a single tile.

        Parameters
        ----------
        mc : ModelContainer
            A ModelContainer containing the input images for a single tile.

        skystat : str
            The sky statistics method to use for skymatching.  Must be 'mean',
            'median', or 'mode'.  See notes.

        Returns
        -------
        ImageModel
            The drizzled 2D image result.

        Example
        -------
        mc = ModelContainer([im1, im2, im3])  # Create a ModelContainer with input images
        skystat = 'median'  # Specify the sky statistics method
        i2d = driz_tile(mc, skystat)
        # `i2d` now contains the drizzled 2D image result.

        Notes
        -----
        - This function configures and runs the JWST calwebb_image3.Image3Pipeline
          on the input ModelContainer for a single tile.  Specifically, it
          runs just the skymatch step and resample step.
        - Within a visit, the skymatch step works decently well, so we can have it
          compute the offsets between exposures in a single visit.
        - While skystat is a parameter, it appears the median works the best
        - This function requires a reference gWCS (to ensure all tiles are drizzled
          to the same grid).  This is stored in 'miri_dummy_wcs.asdf'.  See
          _create_output_grid for details.
        """
        im3 = calwebb_image3.Image3Pipeline()

        # Configure pipeline steps
        im3.tweakreg.skip = True
        im3.outlier_detection.skip = True
        im3.source_catalog.skip = True

        im3.skymatch.skystat = self.skystat
        im3.skymatch.lsigma = self.lsigma
        im3.skymatch.usigma = self.usigma
        im3.skymatch.binwidth = self.binwidth
        im3.skymatch.match_down = self.match_down

        im3.resample.kernel = 'square'
        im3.resample.weight_type = 'exptime'
        im3.resample.save_results = self.save_tiles
        if self.save_tiles:
            im3.resample.output_file = f'{mc.meta.filename}_i2d.fits'

        # Run the skymatch step
        mc = im3.skymatch.run(mc)

        # Configure and run the resample step
        im3.resample.output_wcs = self.wcs_filename
        i2d = im3.resample.run(mc)

        return i2d

    def combine_with_exp(self, tile, tile_bg):
        """
        Adds the tile-to-tile background to the backgrounds computed for each exposure
        within the tile.  This is needed for later redrizzling all the tiles input
        exposures together.

        Parameters
        ----------
        tile : ImageModel
            A JWST Level 3 ImageModel from which to extract the input exposures data.

        tile_bg : float
            The background value to be added to the 'BKGLEVEL' column in the header table.

        Returns
        -------
        astropy.table.Table
            A table containing the combined data with the 'FILENAME' and 'MATCHEDBG' columns.


        Notes
        -----
        - The input `tile` should be a JWST Level 3 ImageModel object.
        - The 'BKGLEVEL' values in the header table of the `tile` object will be
          combined with the `tile_bg` value to create the 'MATCHEDBG' column in the result.
        """
        hdrtab = Table(tile.hdrtab)

        if 'BKGLEVEL' not in hdrtab.colnames:
            hdrtab['BKGLEVEL'] = 0.

        nanmask = np.isnan(hdrtab['BKGLEVEL'])
        hdrtab['BKGLEVEL'][nanmask] = 0.
        hdrtab['MATCHEDBG'] = hdrtab['BKGLEVEL'] + tile_bg
        return hdrtab['FILENAME', 'MATCHEDBG']

    def group_exps_by_visit(self, input):
        """
        Group JWST ModelContainer objects by observation visit.

        Parameters
        ----------
        big_container : ModelContainer
            A JWST ModelContainer object containing data from multiple observation visits.

        Returns
        -------
        list of ModelContainer
            A list of ModelContainer objects, each containing data from a single observation visit.

        Example
        -------
        big_container = ModelContainer([im1, im2, im3, im4])  # Create a ModelContainer with images from multiple visits
        visit_containers = group_exps(big_container)
        # `visit_containers` now contains ModelContainer objects grouped by observation visit.

        Notes
        -----
        - This function groups images from the input `big_container` by their observation visit.
        - Each resulting ModelContainer contains data from a single observation visit.
        """
        all_visits = [im.meta.observation.visit_id for im in input]
        all_dets = [im.meta.instrument.detector for im in input]
        visits = sorted(list(set(all_visits)))
        dets = sorted(list(set(all_dets)))
        mclist = []

        for visit in visits:
            # Create a list of images for the current visit
            visit_images = [im for im in input if im.meta.observation.visit_id == visit]
            for det in dets:
                det_images = [im for im in visit_images if im.meta.instrument.detector == det]
                # Create a ModelContainer for the current visit/detector and append it to mclist
                if det_images:
                    mclist.append(ModelContainer(det_images))

        return mclist

    def clipmed(self, vals):
        """Computes sigma clipped median of values"""
        vals = vals[(np.isfinite(vals))&(vals!= 0.)]
        if len(vals)==0:
            raise IndexError('No values to calculate clipped median.  Ensure arrays overlap')
        clip = sigmaclip(vals.ravel(), self.lsigma, self.usigma)[0]
        return np.median(clip)

    def pair_diffs(self, d1, d2):
        """
        Computes the sigma clipped median difference between two overlapping arrays.

        Parameters
        ----------
        d1 : numpy.ndarray
            The first input array for computing differences.

        d2 : numpy.ndarray
            The second input array for computing differences.

        Returns
        -------
        float
            The clipped median difference between the non-zero, finite elements of `d1` and `d2`.

        Notes
        -----
        - The function computes the clipped median difference between the elements of `d1` and `d2`
          that are both non-zero and finite.
        - Elements with the value 0.0 are excluded from the computation.
        """
        mask = np.isfinite(d1) & np.isfinite(d2) & (d1 != 0.0) & (d2 != 0.0)
        return self._skystat(d1[mask] - d2[mask])

    def _calc_overlap(self, d1, d2):
        """Calculates number of overlapping nonzero, finite pixels in two arrays"""
        mask = np.isfinite(d1) & np.isfinite(d2) & (d1!=0.) & (d2!=0.)
        return np.sum(mask)

    def _create_matrix(self, tiles):
        """Creates the system of equations for the background offsets"""
        pdiffs = []
        eqns = []
        true_overlap = []

        nt = len(tiles)
        print('calculating overlaps')
        for i in range(nt):
            for j in range(i+1, nt):
                # Even though _skystat calculates overlap it raises an error
                # if pixels dont overlap, so need to calculate manually
                overlap = self._calc_overlap(tiles[i].data, tiles[j].data)
                if overlap > 0:
                    skystat_result = self.pair_diffs(tiles[i].data, tiles[j].data)
                    pdiffs.append(skystat_result[0])
                    true_overlap.append(skystat_result[1])
                    row = np.zeros(nt).astype(np.float64)
                    row[i] = 1.
                    row[j] = -1.
                    eqns.append(row)
        return np.array(eqns), np.array(pdiffs), np.array(true_overlap)

    def _set_sky_levels(self, img, bg_dict):
        """Puts the new sky value in the metadata, subtracts sky if set"""

        for dm in img:
            dmname = dm.meta.filename
            dm.meta.background.method = str(self.skymethod)

            sky_level = bg_dict[dmname]
            dm.meta.background.level = sky_level

            if self.subtract:
                dm.data -= sky_level
            dm.meta.background.subtracted = self.subtract

            dm.meta.cal_step.skymatch = "COMPLETE"
