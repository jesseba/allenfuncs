#!/usr/bin/env python
# coding: utf-8

# In[112]:


import os
import os.path as pth
import matplotlib.pyplot as plt
import IPython.display as Disp
from ipywidgets import widgets
import numpy as np
from skimage.transform import PiecewiseAffineTransform, warp
from skimage.transform import PolynomialTransform, warp
from skimage import data
from skimage import io
from scipy.interpolate import interp2d
from PIL import Image as im
from matplotlib import cm
from skimage.transform import rotate
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import pandas as pd
import numpy as np
from matplotlib.colors import ListedColormap
import seaborn as sns

import holoviews as hv
from holoviews import opts
hv.extension('bokeh')


# In[113]:


from allensdk.core.reference_space_cache import ReferenceSpaceCache
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache


# In[114]:


# anatomy
inj = 'retina'
proj = 'Dorsal part of the lateral geniculate complex'

proj_naming = 'dLGN'

#mouse strains used
strain1 = 'Slc17a6-IRES-Cre'
strain2 = 'Foxp2-IRES-Cre'
strain3 = 'Kcng4-Cre' 

#rotation angle
angle = (21.9835,0,-7.2348)


# In[17]:


# get struct IDs for inj/proj, make masks for projection, find experiment ids based on strain and inj/proj, 
# go through IDs, incorporate angle rotation and make zstack ... of just the mask? find where mask is not 0... 
# using zstack, put slice_idx to make some figures


# In[18]:


reference_space_key = 'annotation/ccf_2017'
resolution = 25
rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
#ID 1 = adilt mouse structure graph
tree = rspc.get_structure_tree(structure_graph_id=1)

annotation, meta = rspc.get_annotation_volume()
# The file should appear in the reference space key directory
os.listdir(reference_space_key)

rsp = rspc.get_reference_space()
mcc = MouseConnectivityCache()
all_experiments = mcc.get_experiments(dataframe=True)
structure_tree = mcc.get_structure_tree()


mako = cm.get_cmap('mako',256)
rocket = cm.get_cmap('rocket',256)
viridis = cm.get_cmap('viridis',256)


# In[19]:


def find_ids(inj, proj, strain):
    """""
    inj = injection site exactly how it's written for the Allen ex: 'retina'
    proj = projection area of interest, ex: 'Dorsal part of the lateral geniculate complex'
    strain = mouse lines, 'False' if 'wild' otherwise ex: 'FoxP2-IRES-Cre'
    
    use function to find the IDs associated with the structures from the mouse connectivity cache in order to 
    quickly find experimental IDs relevant to the injection sites and mouse lines of interest
    
    """""
    
    # get tree structures for anatomical regions
    proj_struct =  tree.get_structures_by_name([proj])
    inj_struct = tree.get_structures_by_name([inj])
    
    #get IDs for structures
    struct_id_proj = proj_struct[0]['id'] 
    struct_id_inj = inj_struct[0]['id'] 
    
    #find experiment #s based on strain and proj/inj 
    structures = structure_tree.get_structures_by_name([inj, proj])
    
    if strain == 'wild':
        experiments = mcc.get_experiments(cre=False, injection_structure_ids=[struct_id_inj])
        exp_id = [ e['id'] for e in experiments ]
    else:
        experiments = mcc.get_experiments(cre=[strain], injection_structure_ids=[struct_id_inj])
        exp_id = [ e['id'] for e in experiments ]
    
    output = {'injection_id':struct_id_inj,
              'projection_id':struct_id_proj,
              'experiment_id':exp_id
             }
    
    return output
    


# In[20]:


def mask_on(inj,proj, strain):
    """""
    inj_id = injection ID found through find_ids function
    proj_id = projection ID found through find_ids function
    
    use function to create boolean mask for structures in the Allen and the indices where the structures exist in
    the 3D matrix of the brain (returns 3D matrix or 3 arrays - 0: coronal, 1: transverse, 2: sagittal)
    """""
    
    ids = find_ids(inj,proj,strain)
    proj_id = ids['projection_id']
    inj_id = ids['injection_id']
    
    #make a mask for the structure
    proj_mask = rsp.make_structure_mask([proj_id]) #coronal, transverse, sagittal
    
    #find indices where these masks are active
    proj_idx = np.nonzero(proj_mask) # returns 3 arrays for each dimension
    
    output = {'projection_mask': proj_mask,
              'proj_mask_idx':proj_idx
             }
    return output


# In[21]:


def projden(inj, proj, strain):
    """""
    inj = injection site exactly how it's written for the Allen ex: 'retina'
    proj = projection area of interest, ex: 'Dorsal part of the lateral geniculate complex'
    strain = mouse lines, 'False' if 'wild' otherwise ex: 'FoxP2-IRES-Cre'
    
    Use this function to create 3D projection density matrix
    """""
    
    #get IDs
    ids = find_ids(inj, proj, strain)
    proj_id = ids['projection_id']
    exp_ids = ids['experiment_id']
    
    #get proj masks
    masks = mask_on(inj, proj, strain)
    proj_mask = masks['projection_mask']
    proj_mask_idx = masks['proj_mask_idx']
    
    #get projection densities for all experiment IDs 
    shape = proj_mask.shape + (len(exp_ids),)
    pd = np.empty(shape=shape)
    pd_info = []
    
    for i, x in enumerate(exp_ids):
        pd[:,:,:,i], pd_info_copy = mcc.get_projection_density(x)
        pd_info.append(pd_info_copy)
        
    #norm pd
    pd_norm = np.empty(shape=pd.shape)
    
    for i in range(pd.shape[3]):
        pd_norm[:,:,:,i] = pd[:,:,:,i]/pd[:,:,:,i].max()
    
    #template and annotation
    template, template_info = mcc.get_template_volume()
    annot, annot_info = mcc.get_annotation_volume()
    
    #average the pd across all experiments
    pd_avg = np.mean(pd,axis=3)
    
    #apply masks to pd
    pd_mask = np.empty(shape=shape)
    
    for i, x in enumerate(exp_ids):
        pd_mask[:,:,:,i] = np.multiply(proj_mask,pd[:,:,:,i])
        
    #apply masks to pd_avg
    pdavg_mask = np.multiply(proj_mask,pd_avg)
    
    output = {'pd_info':pd_info,
              'experiment_ids':exp_ids,
              'pd_all_exp' : pd,
              'pd_norm' : pd_norm,
              'pd_all_mask' : pd_mask,
              'pd_avg_exp' : pd_avg,
              'pd_avg_mask' : pdavg_mask,
              'template' : template,
              'annotation' : annot
             }
    
    
    return output


# In[28]:


def rot_projden(inj, proj, strain, angle=(0,0,0)):
    """""
    inj = injection site exactly how it's written for the Allen ex: 'retina'
    proj = projection area of interest, ex: 'Dorsal part of the lateral geniculate complex'
    strain = mouse lines, 'False' if 'wild' otherwise ex: 'FoxP2-IRES-Cre'
    angle = rotation angle if necessary - default is non-rotated (roll, yaw, pitch)
    
    Use this function to rotate the projection density coronally (for window, GRIN and cannula
    implants)
    """""
    
    proj_den = projden(inj,proj,strain)
    pd_info = proj_den['pd_info']
    pd = proj_den['pd_all_exp'] #4D matrix
    pd_norm = proj_den['pd_norm'] #4D matrix normed
    pd_avg = proj_den['pd_avg_exp'] #3D average across experiments
    template = proj_den['template'] 
    annot = proj_den['annotation']
    
    pd_shape = pd.shape #empty 4D matrix
    temp_shape = template.shape

    #apply rotation to normed 4D matrix
    rot_pd_norm = np.empty(shape= pd_shape)
    rot_temp = np.empty(shape= temp_shape)
    
    for i in list(range(pd.shape[0])):
        rot_pd_norm[i,:,:,:] = rotate(pd_norm[i,:,:,:],angle[0])
        rot_temp[i,:,:] = rotate(template[i,:,:],angle[0])
    
    for i in list(range(pd.shape[1])):
        rot_pd_norm[:,i,:,:] = rotate(rot_pd_norm[:,i,:,:],angle[1])
        rot_temp[:,i,:] = rotate(rot_temp[:,i,:],angle[1])
        
    for i in list(range(pd.shape[2])):
        rot_pd_norm[:,:,i,:] = rotate(rot_pd_norm[:,:,i,:],angle[2])
        rot_temp[:,:,i] = rotate(rot_temp[:,:,i],angle[2])
        
    #apply rotation to mask
    mask = mask_on(inj,proj,strain)
    proj_mask = mask['projection_mask']
    
    rot_mask = np.empty(shape=proj_mask.shape)
    
    for i in list(range(proj_mask.shape[0])):
        rot_mask[i,:,:] = rotate(proj_mask[i,:,:],angle[0])
    
    for i in list(range(proj_mask.shape[1])):
        rot_mask[:,i,:] = rotate(rot_mask[:,i,:],angle[1])
        
    for i in list(range(proj_mask.shape[2])):
        rot_mask[:,:,i] = rotate(proj_mask[:,:,i],angle[2])
    
    #rotated mask idx
    rot_mask_idx = np.nonzero(rot_mask)
      
    # average the rotation 4D matrix across experiments
    rot_pd_avg = np.mean(rot_pd_norm,axis=3)
    
    output = {'proj_den' : proj_den,
              'maskon' : mask,
              'rot_pd_norm' : rot_pd_norm,
              'rot_pd_avg' : rot_pd_avg,
              'rot_temp' : rot_temp,
              'rot_mask' : rot_mask,
              'rot_mask_idx': rot_mask_idx
             }
    
    return output


# In[46]:


def flip_rot(inj, proj, strain, angle=(0,0,0)):
    
    rot_projden_out = rot_projden(inj, proj, strain, angle)
                              
    rot_pd_norm = rot_projden_out['rot_pd_norm'] #normalized rotated 4D mat                       
    projden_out = rot_projden_out['proj_den']
    ids = projden_out['experiment_ids']
    pd_norm = projden_out['pd_norm'] #normalized 4D mat
    
    #flip these experiments
    pd_norm_flp = np.empty(shape=pd_norm.shape)
    
    for i in range(pd_norm.shape[3]):
        pd_norm_flp[:,:,:,i] = np.flip(pd_norm[:,:,:,i], axis=2)
    
    ## rotate the other hemisphere

    rot_pd_norm_flp = np.empty(shape= pd_norm_flp.shape)

    for x in list(range(pd_norm_flp.shape[3])):
        
        for i in list(range(pd_norm_flp.shape[0])):
            rot_pd_norm_flp[i,:,:,x] = rotate(pd_norm_flp[i,:,:,x],angle[0])

        for i in list(range(pd_norm_flp.shape[1])):
            rot_pd_norm_flp[:,i,:,x] = rotate(rot_pd_norm_flp[:,i,:,x],angle[1])

        for i in list(range(pd_norm_flp.shape[2])):
            rot_pd_norm_flp[:,:,i,x] = rotate(rot_pd_norm_flp[:,:,i,x],angle[2])
        
    
    output = {'ids' : ids,
              'rot_projden':rot_projden_out,
              'projden' : projden_out,
              'pd_norm' : pd_norm,
              'pd_norm_flp' : pd_norm_flp,
              'rot_pd_norm' : rot_pd_norm,
              'rot_pd_norm_flp' : rot_pd_norm_flp
             }
    
    return output

