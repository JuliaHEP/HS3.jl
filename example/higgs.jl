using HS3
using BenchmarkTools
using DensityInterface

#dict = file_to_dict("/net/e4-nfs-home.e4.physik.tu-dortmund.de/home/rpelkner/Documents/higgs/br_test/llll/125.json")
dict = file_to_dict("/net/e4-nfs-home.e4.physik.tu-dortmund.de/home/rpelkner/Documents/higgs/br_test/llll/126.json")


specs = HS3.generate_specs(dict)

