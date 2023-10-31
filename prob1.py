import json
import numpy as np
from scipy.spatial.transform import Rotation

def transform_coordinates(landmarks, theta, t):
    """Apply the rotation and translation to the given landmarks."""
    r = Rotation.from_euler('xyz', theta, degrees=True)  # Create rotation object from Euler angles
    matrix = r.as_matrix()  # Convert the rotation to a matrix
    transformed = np.dot(landmarks, matrix.T) + t
    return transformed

def main():
    # load json
    with open('Project4.json') as f:
        data = json.load(f)

    Proj4landmarks = np.array(data['Proj4landmarks'])
    
    # Transformation parameters
    theta1 = [10, -15, 5]
    theta2 = [5, 10, 2]
    theta3 = [-5, -5, -7]
    t1 = [2, -1, 4]
    t2 = [-1, 4, 0]
    t3 = [-0.5, -3, -4]

    # Applying transformations
    Proj4landmarks = transform_coordinates(Proj4landmarks, theta1, t1)
    Proj4landmarks = transform_coordinates(Proj4landmarks, theta2, t2)
    Proj4landmarks = transform_coordinates(Proj4landmarks, theta3, t3)

    # numpy round all values to 3 decimal places
    Proj4landmarks = np.round(Proj4landmarks, 0)

    # cast to int
    Proj4landmarks = Proj4landmarks.astype(int)

    print(Proj4landmarks)

if __name__ == "__main__":
    main()
