import subprocess
import os

print 'Suffix: '
suffix = raw_input(">>")

directory_path  = '/KM/disopt/archive00/sar/'
executables     = ['../measure_lemon_cost_scaling']
sub             = 'generated_grids_dimacs/'
timeout         = 10000
number_of_runs  = 5
# suffix          = 'a'
log_name        = '../../data/results/benchmarks/cosgeneratedgrids_' + suffix + '.txt'
lis = [graph for graph in os.listdir(directory_path + sub) if '_256_' + suffix + '_' in graph and graph.endswith('.min')]
lis = sorted(lis)

for executable in executables:
  print executable
  with open(log_name,'a') as f:
    print >> f, 'Running ' + executable + ' ' + str(number_of_runs) + ' times on graphs in ' + sub + ' with timeout ' + str(timeout) + '\n'
  for graph in lis:
    with open(log_name,'a') as f:
      print >> f, 'new graph: ', graph, '\n'
    for i in range(number_of_runs):
      cmd = executable + ' ' + directory_path + sub + graph + ' ' + str(timeout) + ' 2>&1'
      print cmd
      output = subprocess.check_output(cmd +  '; exit 0', stderr=subprocess.STDOUT, shell=True)
      with open(log_name,'a') as f:
        print >>f, output
      if 'Timeout!' in output:
        break
