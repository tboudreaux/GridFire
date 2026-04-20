import pynucastro as pyna

def main():
    lib = pyna.ReacLibLibrary()
    
    filtered_rates = []
    
    for rate in lib.get_rates():
        if all(nuc.Z <= 26 for nuc in rate.reactants):
            filtered_rates.append(rate)
            
    out_file = "reaclib_pynucastro_latest_Z26.dat"
    
    with open(out_file, "w") as f:
        for rate in filtered_rates:
            rate.write_to_file(f)
            

if __name__ == "__main__":
    main()
