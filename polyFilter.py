import sys
from polyLine import polyLine




if __name__ == "__main__":
    polyFile = sys.argv[1]

    with open(polyFile) as file:
        for line in file:
            pLine = polyLine(line)
            sys.exit()




