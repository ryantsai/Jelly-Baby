# Bed and sleeping

Both blanket contact passes use the existing soft-body WASM module when
available, preserving the same triangle traversal and projection order.
See [Native collision kernels](native-collision.md) for shared-memory ownership
and exact JavaScript comparison checks.

The bed occupies the foreground below the trampoline: world `(.165, .205)` in
X/Z. Its long axis follows Z, with the pillow at negative Z. Rounded oak,
an inset headboard, brass foot collars, mattress, piped pillow, and sage
coverlet share the existing facility shadow paths. Both maps share the dynamic
fabric geometry, supporting jelly-to-blanket and blanket-to-jelly occlusion.

Boarding requires grounded, ungrabbed proximity. A slightly inclined supine
frame supports the head over the pillow and leaves the face above the coverlet.
Damped per-node support and gravity compensation retain the ordinary FEM
deformation and volume solve. Small chest support motion supplies breathing.
Getting up reverses the frame onto the clear left approach and clears velocity.

The blanket is a 57 by 49 height-field membrane. Gravity, neighboring tension,
and damping evolve vertical motion at the fixed physics rate. Textile X/Z
coordinates remain fixed, excluding tangling and lateral sliding. The mattress
supplies lower contact heights, including hanging side hems. The bound optical
surface proxy is reconstructed after physics and its triangles rasterized into
vertical contact heights with a 0.65 mm fabric allowance. Contact footprints
are not expanded: spreading neighboring maxima visibly inflated the opening
around the body. A light damped conforming preload removes excess lift while
retaining breathing response. The finer grid follows the torso and arms without
broad tension bridges.

Proxy contact alone is insufficient at steep sides: a fabric triangle can cut
through a curved arm even when its grid vertices are clear. After the visible
skin is updated, a separate render contact pass clips its exact triangles
against the fabric triangles in X/Z. Each intersection polygon vertex imposes
a 0.45 mm clearance constraint. Since both surfaces are piecewise linear,
these constraints also cover their overlapping interiors. Corrections only
raise low fabric corners, preserving earlier contacts without inflating an
already-clear crown. This pass uses reusable buffers once per rendered frame,
not per fixed physics step. The corrected mesh is shared with both shadow paths.
Corrections do not feed displacement or energy into the membrane simulation.
Contact corrections can leave grid-aligned ridges on the hanging skirt. The
top should fit the jelly, but the hidden underside does not need to trace every
body contour. A bounded tension relaxation spans roughly two centimetres of
fabric (up to 128 symmetric passes with an early settling threshold), letting
the sides bridge away from the body into a broad drape. The earlier eight-pass
local smoothing retained too much of the tight ridge. Relaxation raises valleys
above the contact-safe envelope, so clearance survives; convex supported tops,
the fitted opening, and outer hem retain their support. Cloth shading normals use symmetric grid tangents
instead of triangle-area weighting, avoiding diagonal lighting bias where
small crown triangles meet tall side triangles. The fairing affects actual
geometry and therefore its shadows, rather than masking folds in the material.
Regression checks also bound discrete crease energy and the steepest interior
grid drop, and verify that relaxation never lowers the contact envelope.

The blanket opts into a curved receiver depth map in `SurfaceShadows`. Wide
shadow-filter taps cannot extrapolate the current triangle's tangent plane
across a curved skirt: that comparison makes neighboring parts of the same
surface appear to be blockers. The additional depth map stores the actual
blanket geometry in the same light projection. Each tap compares against that
receiver depth; only the subtexel center correction uses the tangent plane.
The center anchor takes the conservative envelope of four surrounding texels,
avoiding nearest-texel jumps where neighboring texels belong to different
triangles. Finite-resolution grazing triangles can still differ from continuous
visibility; the regression bounds this subpixel disagreement below 0.1% rather
than asserting exact visibility everywhere.
The center's separation from the nearest blanket layer is retained, preserving
occlusion by overlapping folds. Other facilities and the jelly still cast onto
the blanket, and the blanket remains in both existing caster paths. The map
tracks geometry, transforms, visibility, and lighting, and is disposed with the
shadow system. `npm run test:blanket-shadows` rasterizes the actual drape on the
CPU to reproduce the old false shadows and checks the new comparison, retained
blockers, caching, and lighting invalidation.

Render geometry versions track the corrected positions, including changes
caused by skin deformation. Regression checks sample all visible skin vertices,
edge midpoints, and triangle centers throughout breathing, while retaining the
close-fit check at the opening. Coupling is one-way:
the lightweight blanket does not load the jelly. This bounded approximation
does not attempt arbitrary cloth folds. Geometry versions change only when
the membrane moves, retaining shadow caching after settling.

Sleep blends curved closed eyelids, lowered brows, a small mouth, and a fading
tongue through the original skin attachment. A nose bubble has a narrowed
attachment pole, damped inflation, and velocity-dependent shape lag. Its
transparent physical material supplies wet highlights and subtle thin-film
color without an extra scene refraction pass. Wake reverses the blend and
deflates the bubble. Reset clears the bubble immediately.

The facility overrides the shared prompts: desktop uses `Press E to Go to Bed`
and `Press E to Get Up`; touch uses `Go to Bed` and `Get up`. Ownership, camera,
input suppression, and teardown use the existing facility manager.

`npm run test:physics` includes the modular bed regression after the core physics
checks. Its facial-clearance check covers skin ink separately from the bubble,
whose geometry is positioned by an object transform.

`npm run test:bed` also runs independently and checks boarding, ten seconds of occupied FEM simulation,
volume, blanket lift and settling, whole-surface clearance, sampled face updates, wake/reset, and
shadow flags and bounds. Numerical checks do not establish visual quality;
browser inspection remains a separate manual check.
