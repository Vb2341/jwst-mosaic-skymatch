jwst-mosaic-skymatch
======================

This package contains tools for improving the background matching for JWST
pipeline produced mosaics.  Currently, it contains a single tool to optimize
sky matching for MIRI image data, as the current default pipeline often
produces level 3 products with discontinuities/artificial gradients in the
background.

Specifically, this package contains a replacement step for the `SkyMatchStep`
of the official JWST pipeline.  The replacing step, `PixelSkyMatchStep`, first
groups all of the exposures by their visit ID's (as images in the same visits
tend to overlap).  For each group of exposures, the exposures are drizzled
(resampled) together but projectedto the same pixel grid for all the groups.
This allows for offsets in the background to be calculated per pixel, leading
to better measurements of background offsets.  These calculated offsets are
then placed back into the exposure metadata, and so the data can be processed
through the remaining pipeline steps (outlier detection and resample) as
usual.

Caveats:
*This algorithm requires the exposures to be aligned to well within a pixel, or
the calculated background offsets will be incorrect.
*This replacement step does not yet support all of the different options as
seen in the default pipeline step.  It currently only replaces the case in which
the sky is desired to be matched (i.e. `skymethod='match'`), though the
global+match method will likely be implemented in the future.
*Only a few dither strategies have been tested with this algorithm, and so other
dither combinations may not see much improvement (TBD).
*This is really only meant for MIRI imaging data.  Other instruments may not
work.

License
-------

See `LICENSE.rst` for more information.


Contributing
------------

We love contributions! `jwst-mosaic-skymatch` is open source,
built on open source, and we'd love to have you hang out in our community.

**Imposter syndrome disclaimer**: We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

*Note:* This disclaimer was originally written by
`Adrienne Lowe <https://github.com/adriennefriend>`_ for a
`PyCon talk <https://www.youtube.com/watch?v=6Uj746j9Heo>`_, and was adapted by
packagename based on its use in the README file for the
`MetPy project <https://github.com/Unidata/MetPy>`_.
