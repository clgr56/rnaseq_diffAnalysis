import glob
import subprocess

def get_names():
    names=[]
    for fname in glob.glob("data/Berlin_AKI_tubuloid_bulk_RNA_seq/*.gz"):
        names.append(fname)
    names_out=[]
    for n in names:
        n_out = n + 'rrna.fq'
        names_out.append(n_out)
    return names, names_out

def rrna(input,output):
    cmd = ['ribodetector_cpu', '-t', '20', '-l', '150', '-i', '-e', 'rrna', '--chunk_size', '256', '-o']
    cmd.insert(6,input)
    cmd.extend(output)
    subprocess.run(cmd,check=True)

if __name__=="__main__":
    n_in, n_out = get_names()
    rrna(n_in,n_out) 
