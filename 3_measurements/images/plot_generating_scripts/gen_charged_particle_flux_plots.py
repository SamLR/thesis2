"""
A script that generates the 1 and 2D charged particle flux plots
"""

from ValueWithError import ValueWithError

count_time_s = 50

class DataPoint(object):
  """docstring for dataPoint"""
  def __init__(self, x_pos, y_pos, hit_count, intensity, notes=None, count_time=count_time_s):
    super(dataPoint, self).__init__()
    self.x_pos = x_pos
    self.y_pos = y_pos
    self.hit_count  = ValueWithError(hit_count)  
    self.intensity  = ValueWithError(intensity)
    self.count_time = ValueWithError(count_time)
    self.notes = notes
    
  def get_normalised_rate():
    
    
count_bias     = (30.0 + 83.0 + 34.0)/3.0
intensity_bias = (243.0 + 381.0 + 398.0 + 404.0)/4.0

data_1D = (DataPoint(6, -15,  446302, 1212),
           DataPoint(6, -15,  502596, 1208),
           DataPoint(6,  -5, 1438685, 1286),
           DataPoint(6,   0, 2151736, 1255),
           DataPoint(6,   0, 1307015, 1398, "MPPC3: dead, MPPC2: damaged"),
           DataPoint(6,   5, 1663702, 1298, "MPPC3: dead, MPPC2: damaged"),
           DataPoint(6,   5, 1420836, 1142, "MPPC3: dead, MPPC2: damaged"),
           DataPoint(6,  15, 1080170, 1336, "MPPC3: dead, MPPC2: damaged"),
           DataPoint(6,  15, 1185051, 1371, "MPPC3: dead, MPPC2: damaged"))


def main():
  # TODO Finish this
  

if __name__=="__main__":
  main()