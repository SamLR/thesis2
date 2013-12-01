"""
Process the efficiency file, then calculate efficiencies. 
"""

from datetime import datetime
from ValueWithError import ValueWithError as Value
from math import pi

def calculateRunTime(raw_data, time_format='%y-%m-%d %H:%M'):
    """
    Pop the start/stop datetimes, convert them and calculate the 
    difference, then add it to the data dictionary.
    """
    start = raw_data.pop('start_date') + " " + raw_data.pop('start_time')
    stop = raw_data.pop('stop_date') + " " + raw_data.pop('stop_time')
    start = datetime.strptime(start, time_format)
    stop = datetime.strptime(stop, time_format)
    # run_time is now a timedelta, error on measurement is ~0.5 minute
    raw_data['run_time'] = Value((stop-start).total_seconds(), 30)

def readFile(file_name):
    data_format = ('start_date', 'start_time', 'stop_date', 'stop_time',
                    'thres', 'config', 'mppc_count', 'pmt_count')
    res = {}
    with open(file_name, 'r') as in_file:
        for line in in_file:
            # skip comments
            if line.startswith('#'): continue
            # split up the table into colums
            tokens = [t.strip() for t in line.split("|")]
            # create our dictionary 
            data = dict(zip(data_format, tokens))
            calculateRunTime(data)
            data.pop('thres') # pop this as we don't need it
            data['mppc_count'] = float(data['mppc_count'])
            data['pmt_count'] = float(data['pmt_count'])
            data['mppc_rate'] = Value(data['mppc_count'])/data['run_time']
            data['pmt_rate'] = Value(data['pmt_count'])/data['run_time']
            data['eff'] = data['mppc_rate']/data['pmt_rate']
            key = data.pop('config')
            res[key] = data
    return res

def getIndividualEff(total_rate, data):
    return total_rate/(data['mppc_rate']/data['pmt_rate'])

def processData(data):
    tmp = {k:v for (k,v) in data.items()} # lazy copy
    all = tmp.pop('123')
    res = {'raw_eff': (all['mppc_rate']/all['pmt_rate']), 
           'adjusted_eff': (all['mppc_rate']/all['pmt_rate']) * (40*40/(pi*3.5**2))}
    # only 12, 13 and 23 remain now
    for config, d in tmp.items():
        # figure out which MPPC number is _not_ in the key (i.e. key = '12' mppc_id = '3')
        mppc_id = filter(lambda x: x not in config, ('1', '2', '3'))[0]
        res[mppc_id] = all['eff']/d['eff']
        
    return res

def prettyPrint(efficiencies):
    pass

def main():
    data = readFile('music_2_efficiency.txt')
    for i in data: 
        print i, data[i]['run_time'], data[i]['mppc_count'], data[i]['pmt_count'], 
        print data[i]['mppc_count']/data[i]['pmt_count'], data[i]['mppc_count']/data[i]['pmt_count']* (40*40/(pi*3.5**2))
    print
    efficiencies = processData(data)
    for e in efficiencies: print e, efficiencies[e]
    prettyPrint(efficiencies)

if __name__ == '__main__':
    main()