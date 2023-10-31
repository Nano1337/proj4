from ctypes import*
import numpy as np

class ndi:
    def __init__(self, ConfigFile):
        self.mydll = cdll.LoadLibrary('C:\\Users\\noblejh\\Box Sync\\Code\\NDI\\x64\\Release\\ndi.dll')
        self.mydll.Init.restype = c_void_p
        self.h = c_void_p(self.mydll.Init(create_string_buffer(ConfigFile)))

    def Start(self):
        self.mydll.Start(self.h)

    def Stop(self):
        self.mydll.Stop(self.h)

    def GetNumTools(self):
        n = self.mydll.GetNumTools(self.h)
        return n.value

    def GetTools(self,n):
        ToolMats_c = c_double(np.zeros([4,4,n]))
        self.mydll.GetTools(self.h,ToolMats_c)
        return ToolMats_c.value

    def __del__(self):
        self.mydll.Delete(self.h)
        self.h = 0