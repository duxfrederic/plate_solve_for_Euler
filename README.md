# a plate solver for Coralie's bigeye acquisition camera

A package containing the following logic:
- acquisition image contains a lot of stars: plate solve acquisition image to find target
- acquisition image contains a few (< 6) stars with similar fluxes: cannot automatically handle this case for now, fall back to manual
- acquisition image contains a few (< 6) stars with one much brighter than the others: given that expolanet stars targets are typically bright, go for it. The PI should specify a "manual" acquisition if they know their star is next to a much brighter one.
- acquisition image contains only one bright star: go for that one.


## TODO
- for now we are using the API at nova.astrometry.net for plate solving. This does not work in production, so install astrometry.net on the glsguiding machine. (use the 4100 and 5200 index files)
- write the plate solve function using the local installation
- convert the J2000 coordinates offsets into new alt az coordinates the telescope should slew to. If we insist on sending pixel offsets, convert with invert gnomonic projection.
- 