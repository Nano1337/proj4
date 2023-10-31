import numpy as np
import json

from DisplayVolume import DisplayVolume, imagevolume

def imageProjection(source, targetdim, sourcevoxsz=[], targetvoxsz=[], T=[], D=[]):
    # Initialize the projected image with zeros
    projected = np.zeros(targetdim)
    
    # Check if we are using a transformation matrix
    if len(T) != 0:
        # Go through each voxel in the source image
        for i in range(source.shape[0]):
            for j in range(source.shape[1]):
                for k in range(source.shape[2]):
                    # Convert voxel to mm
                    point = np.array([i * sourcevoxsz[0], j * sourcevoxsz[1], k * sourcevoxsz[2], 1])
                    # Apply transformation
                    transformed_point = np.matmul(T, point)
                    # Convert back to voxel coordinates in the target image
                    x = int(transformed_point[0] / targetvoxsz[0])
                    y = int(transformed_point[1] / targetvoxsz[1])
                    z = int(transformed_point[2] / targetvoxsz[2])
                    # Check if the transformed point is within the target dimensions
                    if 0 <= x < targetdim[0] and 0 <= y < targetdim[1] and 0 <= z < targetdim[2]:
                        projected[x, y, z] = source[i, j, k]
                        
    # If deformation field is provided
    elif len(D) != 0:
        # Go through each voxel in the source image
        for i in range(source.shape[0]):
            for j in range(source.shape[1]):
                for k in range(source.shape[2]):
                    # Fetch the deformation from the deformation field
                    dx, dy, dz = D[:, i, j, k]
                    # Apply deformation
                    x = i + dx
                    y = j + dy
                    z = k + dz
                    # Check if the deformed point is within the target dimensions
                    if 0 <= x < targetdim[0] and 0 <= y < targetdim[1] and 0 <= z < targetdim[2]:
                        projected[int(x), int(y), int(z)] = source[i, j, k]
    
    return projected

if __name__ == "__main__": 

    # load json
    with open('Project4.json') as f:
        data = json.load(f)

    ct = np.array(data['ct']['data'], dtype=np.int16)
    t1 = np.array(data['t1']['data'], dtype=np.int16)
    CT2T1 = np.array(data['CT2T1'], dtype=np.float64)

    # project ct onto t1 with transformation matrix CT2T1
    projected_t1 = imageProjection(ct, t1.shape, sourcevoxsz=data['ct']['voxsz'], targetvoxsz=data['t1']['voxsz'], T=CT2T1)

    print(projected_t1.shape)

    # Visualizing the original and thresholded images
    d = DisplayVolume()
    d.SetImage(t1, data['t1']['voxsz'])
    # d.SetImage(projected_t1, data['t1']['voxsz'])
    # d.AddMask(ct, color=[1, 0, 0], opacity=0.5, label='Original CT Image')
    # d.AddMask(projected_t1, color=[0, 1, 0], opacity=0.5, label='Projected CT onto T1')
    d.Display()


    # d.SetImage(mr, voxsz)
    # d.AddMask(lvr, color=[1, 0, 0], opacity=0.5, label='Ground Truth Liver Segmentation')
    # d.AddMask(segmented, color=[0, 1, 0], opacity=0.5, label='Otsu Segmentation')
    # d.Display()