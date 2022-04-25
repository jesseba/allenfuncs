# allenfuncs
Functions created to use Allen's Mouse Connectivity Cache easy to use. Functions allow you to set injection and projection sites for specific mouse strains and look at subsequent brain regions of interest. 

I specifically pulled the ccf_2017 from Allen for the reference space key and the Mouse Connectivity Cache. In my code, I use retina for the injection site, dorsal part of the lateral geniculate nucleus as the projection site, and I have three mouse strains based on our project direction (Slc17a6-IRES-Cre, FoxP2-IRES-Cre, Kcng4-Cre). I decided to name files using our acronym for dorsal part of the lateral geniculate nucleus, dLGN, as opposed to the Allen's (GENd). That variable isn't necessary for any of the code committed. 

**Functions:**
_* find_ids:_
Input for this function is the injection area, projection area, and strain. The function pulls all experiments that are narrowed down by the params passed through the function. The output is a dict of IDs for the injection site, projection site as well as relevant experiments. These will come in handy for the next few functions. 

_* mask_on:_
Input for this function is the same as the above (injection area, projection area, and strain). The output is a dict containing the 3D matrix of the boolean mask for the projection area of the brain. This allows you to apply it to the projection matrix (as is done in a later function) to narrow the projection down to that particular brain region and removes the other projections found throughout the brain. It also outputs a 3D matrix of indices that the mask exists to make it easier to index into the projection matrix to find the relevant slices. 

_* projden:_
Input for this function is the same as the above (injection area, projection area, and strain). The output for this function is a dict that contains the projection density info, the experiment IDs, a 4D matrix of the projection densities (3D) for all experiment IDs (4th dim), projection density normalized, the projection density with the mask applied, the average of the projection density across all experiments, the averaged projection density with the mask applied, the image registration template, and the annotation template from Allen. 

_* rot_projden:_
Input for this function is the same as the above (injection area, projection area, and strain) but it also includes an angle (roll, yaw, pitch). This is for people who are rotating the mouse's head in a particular way before conducting their experiment. Typically histology done on brains is sliced in a conventional way with no intentional rotations. This allows us to rotate this 3D matrix on whichever axis to allow us to (for example) recreate the a field of view for those imaging through GRIN lenses or windows or cannulas implanted at an angle. This can also be used for patch physiology, etc. The output for this is a dict that gives you the non-rotated projden function's output, non-rotated mask, rotated projection density normalized, rotated projection density averaged across experiments, rotated image registration template, rotated mask, and the rotated mask indices. 

_* flip_rot:_
Input for this function is the same as the above (injection area, projection area, strain, and angle). The idea for this function may be a bit niche for our project's needs. We were seeing the ipsilateral patch in the field of view from our 2P imaging data. So to find the same slice in the Allen dataset, we decided to overlay both the contralateral and ipsilateral projections of vglut2 to dLGN (diff colors). But for this to work, we needed to rotate the ipsilateral side to the same angle as the contralateral side before the overlay. So that's what this function does - it's output is a dict with the nonrotated projectin density, rotated projection density, and rotated flipped version. We used these matrices to add our own cmaps to them and overlay the images using PIL to alpha blend the overlays.


Let me know if you need any clarifications or have any suggestions!
