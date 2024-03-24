import argparse



def main(*args, **kwargs):
    print(*args)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description= "This program generates vector descriptors for input pdb files.")
    parser.add_argument("-i", "--input", nargs=1, help="Input .pdb file")
    parser.add_argument("-o", "--output", nargs=1,help="Output.npy file")
    args = parser.parse_args()
    main(args.name)
