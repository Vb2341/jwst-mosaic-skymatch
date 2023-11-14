#! /usr/bin/env python
"""
A small trick to make calculate separate sky values for each NIRCam Chip.



:Authors: Varun Bajaj


"""

from jwst.datamodels import ModelContainer
from jwst.skymatch import SkyMatchStep

class UngroupedSkyMatchStep(SkyMatchStep):
    reference_file_types = []
    spec = SkyMatchStep.spec
    def process(self, input):
        """
        This step just runs the SkyMatchStep, but allows every NIRCAM chip to have its own sky value.  See jwst.skymatch.SkyMatchStep for details.

        By default, if all NIRCAM detectors are used simultaneously, all 4 shortwave exposures
        from each module get treated as one image for the purposes of sky matching.  That
        means only a single sky value gets computed for all 4 chips, which only works if
        the 4 chips have the same exact background beforehand.   In practice they do not.
        Overriding the ``exposure_number`` attribute allows each detector to be treated on its own.
        """
        # Make sure it's a model container
        img = ModelContainer(input)
        orig_num = [m.meta.observation.exposure_number for m in img]
        for i, m in enumerate(img):
            setattr(m.meta.observation, 'exposure_number', str(i))

        super().process(input)

        for i, m in enumerate(img):
            setattr(m.meta.observation, 'exposure_number', orig_num[i])

        return img
