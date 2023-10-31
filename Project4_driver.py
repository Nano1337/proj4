import json
import numpy as np
from DisplayVolume import *
import time
from Points_to_Image import *
from mayavi import mlab
import copy

class DisplayReg(DisplayVolume):
    def __init__(self):
        super().__init__()
        self.renwin2 = vtk.vtkRenderWindow()
        self.ren2 = vtk.vtkRenderer()
        self.renwin2.AddRenderer(self.ren2)
        self.renwin2.Render()

    def Update(self,direction=-1,slc=-1,level=float("nan"),contrast=float("nan"),resize=-1):
        super().Update(direction,slc,level,contrast,resize)
        if direction == 3:
            self.ren.ResetCameraClippingRange()
            self.renwin.Render()
            self.inter.ProcessEvents()
            self.UpdateCam2()

    def UpdateCam2(self):
        self.renwin2.MakeCurrent()
        self.ren2.RemoveAllViewProps()
        res = self.ren.GetActors()
        res.InitTraversal()
        a = res.GetNextActor()
        while a is not None:
            self.ren2.AddActor(a)
            a = res.GetNextActor()

        cam = self.ren.GetActiveCamera()
        cam2 = vtk.vtkCamera()
        cam2.DeepCopy(cam)
        cam2.Azimuth(90)
        self.ren2.SetActiveCamera(cam2)

        self.ren2.ResetCameraClippingRange()
        self.renwin2.Render()

    def Display(self, blocking=True):
        self.Update()
        if blocking:
            while (self.quit == False):
                if usemayavi:
                    self.fig.canvas.draw_idle()
                    self.fig.canvas.start_event_loop(0.3)
                else:
                    self.ren.ResetCameraClippingRange()
                    self.renwin.Render()
                    self.inter.ProcessEvents()
                    self.UpdateCam2()
    def Close(self):
        super().Close()
        del self.renwin2,self.ren2


def project4_driver(trk, probepnt):
    with open('Project4_8c.json') as f:
        dt = json.load(f)
    mr = np.array(dt['mr']['data'], dtype=np.int16)
    # lvr = np.array(dt['lvr']['data'], dtype=np.int16)
    Tt1_f2 = np.array(dt['Tt1_f2'], dtype=np.float64)
    frame_V = np.array(dt['frame_V'], dtype=np.float64)
    marker_V = np.array(dt['marker_V'], dtype=np.float64)
    pointer_V = np.array(dt['pointer_V'], dtype=np.float64)
    frame_F = np.array(dt['frame_F'], dtype=np.longlong)
    marker_F = np.array(dt['marker_F'], dtype=np.longlong)
    pointer_F = np.array(dt['pointer_F'], dtype=np.longlong)
    targets = np.array(dt['targets'], dtype=np.float64)

    d = DisplayReg()
    d.SetImage(mr,np.array(dt['mr']['voxsz']))

    # set camera position?


    frame = surface()
    frame.verts = frame_V
    frame.faces = frame_F
    obj = object(0, frame, color=[0.2, 0.2, 0.2], opacity=1)
    d.objs.append(obj)

    pointer = surface()
    pointer.verts = pointer_V
    pointer.faces = pointer_F
    obj = object(0, pointer, color=[0.7, 0.7, 0.7], opacity=1)
    d.objs.append(obj)

    markers = surface()
    markers.verts = marker_V
    markers.faces = marker_F
    obj = object(0, markers, color=[.95, .95, .95], opacity=1)
    d.objs.append(obj)



    n = trk.GetNumTools()

    target_frame = copy.copy(frame)
    target_pointer = copy.copy(pointer)
    target_markers = copy.copy(markers)
    obj = object(0, target_frame, color=[.2, 0, 0], opacity=1)
    d.objs.append(obj)
    obj = object(0, target_pointer, color=[.7, 0, 0], opacity=1)
    d.objs.append(obj)
    obj = object(0, target_markers, color=[.95, 0, 0], opacity=1)
    d.objs.append(obj)
    # d.AddMask(lvr, [0., 1.0, 0.], label='Liver')

    lvrs = surface()
    lvrs.verts = np.array(dt['lvr']['Vertices'])
    lvrs.faces = np.array(dt['lvr']['Faces'])
    obj = object(0, lvrs, color=[1.0, 0.7, 0.7], opacity=1)
    d.objs.append(obj)

    for i in range(np.shape(targets)[2]):
        targethit = 0
        start = time.time()
        while d.quit==0 and targethit==0:
            Tr = trk.GetTools(n)
            testbedPhi = Tr[:, :, 1]
            pointerPhi = Tr[:, :, 2]
            if testbedPhi[0,3] != 0 and pointerPhi[0,3] != 0:
                x = Points_to_Image(Tt1_f2, testbedPhi, pointerPhi, probepnt)

            elif testbedPhi[0,3] == 0:
                print('Demo board blocked from tracker!')
                continue
            elif pointerPhi[0,3] == 0:
                print('Pointer tool not visible by tracker!')
                continue

            d.objs[0].data.verts = Points_to_Image(Tt1_f2, testbedPhi, pointerPhi, frame_V)
            d.objs[1].data.verts = Points_to_Image(Tt1_f2, testbedPhi, pointerPhi, pointer_V)
            d.objs[2].data.verts = Points_to_Image(Tt1_f2, testbedPhi, pointerPhi, marker_V)
            d.objs[3].data.verts = Points_to_Image(Tt1_f2, testbedPhi, testbedPhi @ targets[i, :, :], frame_V)
            d.objs[4].data.verts = Points_to_Image(Tt1_f2, testbedPhi, testbedPhi @ targets[i, :, :], pointer_V)
            d.objs[5].data.verts = Points_to_Image(Tt1_f2, testbedPhi, testbedPhi @ targets[i, :, :], marker_V)

            tpnt = Points_to_Image(Tt1_f2, testbedPhi, testbedPhi @  targets[i, :, :], probepnt)
            targeterror = np.sqrt(np.sum((tpnt[0:3] - x[0:3])**2))
            if targeterror < 10:
                targethit = 1
                if usemayavi:
                    d.objs[3].ms.mlab_source.set(color = (0, 0.2, 0))
                    d.objs[4].ms.mlab_source.set(color = (0, 0.7, 0))
                    d.objs[5].ms.mlab_source.set(color = (0, 0.95, 0))
            else:
                targethit = 0

            print(f'Target error = {targeterror: .1f} mm')
            d.Update(direction=3)
            
        t = time.time() - start
        print(f'You found the target in {t:.1f} seconds with targeterror of' +\
                 f'{targeterror:.1f} mm. Press enter to continue')
        if d.quit:
            break
        input()
        if usemayavi:
            d.objs[3].ms.mlab_source.set(color=(0.2, 0, 0))
            d.objs[4].ms.mlab_source.set(color=(0.7, 0, 0))
            d.objs[5].ms.mlab_source.set(color=(0.95, 0, 0))

