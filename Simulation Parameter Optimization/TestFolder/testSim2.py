import sys
if __name__ == "__main__":
    x = float(sys.argv[1])
    y = float(sys.argv[2])
    z = float(sys.argv[3])
    
    outfile = open("testSimData.txt",'w')
    val = x**4+x**3-4*x**2
    val += y**4+y**3-4*y**2
    val += z**4+z**3-4*z**2
    outfile.write(str(val))




