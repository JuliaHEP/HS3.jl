#Example file for simple HistFactory typed workspaces.

using HS3
using BenchmarkTools, DensityInterface
using JSON3
  
dict = file_to_dict("./example/rf515_hfJSON_3.json")
specs = HS3.generate_specs(dict)
analysis = HS3.make_analyses(specs.analyses[1], specs)
#Make sure to use Floats! Automatic promotion needs to be implemened
logdensityof(analysis[:likelihood], (mu=1., syst1=0., syst2=0.0, syst3 = 0., lumi = 1., gamma_stat_channel1_bin_0 = 0.1, gamma_stat_channel1_bin_1 = 1.))

#Put your specific analysis code here...


### Second workspace

dict = file_to_dict("./example/histfactory_hl.json")
specs = HS3.generate_specs(dict)
analysis = HS3.make_analyses(specs.analyses[1], specs)
#Make sure to use Floats! Automatic promotion needs to be implemened
logdensityof(analysis[:likelihood], (mu=1., alpha_syst1=0., alpha_syst2=0.0, alpha_syst3 = 0., lumi = 1., gamma_stat_channel1_bin_0 = 0.1, gamma_stat_channel1_bin_1 = 1.))

