import subprocess
import os

input_file = '/home/vruizser/PhD/2018-2019/PanCancer_data/vep_output/vep_PanCan_chr_1_1-100000'
out_dir='./scripts/input_pdbmapper/splitted_vcf_db/'
os.makedirs(out_dir, exist_ok=True)


command = "awk -F ' ' '{{for(i=1;i<=NF;i++) {{if ($i ~ /ENSG/){{print i; exit}}}}}}' {} ".format(input_file)


print(command)

# First command
p1 = subprocess.Popen(command, stdout=subprocess.PIPE, shell = True)
# Second command's input linked to the first one's output

# Read from p2 to get the output
out, err = p1.communicate()
#print(p1.stdout)
command2 = "grep -v '##' {} | \
awk -v ci=\"{}\" \
-v od=\"{}\" \
-F ' ' 'NR==1 \
{{h=$0; next}} \
{{f=od$ci\".vep\"}} !($ci in p) \
{{p[$ci]; print h > f}} \
{{print >> f; close(f)}}'".format(input_file, out.decode('utf8'), out_dir) 

#command2 = "grep -v '##' {}".format(input_file) 

print(command2)
p2 = subprocess.Popen(command2, stdout=subprocess.PIPE, shell = True)

# Read from p2 to get the output
out2, err2 = p2.communicate()