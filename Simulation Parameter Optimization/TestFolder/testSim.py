import sys
if __name__ == "__main__":
    x = float(sys.argv[1])
    y = float(sys.argv[2])
    z = float(sys.argv[3])
    
    outfile = open("testSimData.txt",'w')
    val = (x-1.1)**2 + (y+4.232)**2 + (z-0.5)**2
    outfile.write(str(val))




