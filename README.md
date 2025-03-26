# CHARON


# CHARON: Environmental DNA Drift Time and Population Estimation Tool

Environmental DNA or forensic DNA can contain the DNA of multiple human individuals. Sequencing this DNA and aligning it to the human reference will result in alignments in BAM format. 
CHARON: **C**alculating **HA**plotype distribution to **R**esolve **O**riginal **N**umber of individuals and demography takes a BAM file representing whole genome sequencing data from potentially multiple humans and infers:
 
1. **Number of Diploid Individuals**: An estimate of the current population size based on genetic data.
2. **Demography**: Drift time from a particular population and can determine which population has the lowest drift time thus the closest proxy.


## Assumptions

We generally assume
- **equal contribution**: Each individual contrinuted more or less the same amount of DNA to the mix 
- **same population of origin**: All the individuals come from the same population.
  
## Installation

To install the tool, clone this repository and follow the setup instructions provided in the `INSTALL.md` file.

## Usage

1. Prepare your eDNA sample data in the native format.
2. 
3. Run the tool using the command-line interface or graphical user interface.
4. Input your data file when prompted.
5. The tool will process the data and output estimates for drift time and population size.

For detailed instructions and examples, refer to the `USER_GUIDE.md` file.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.


---

Ensure that you replace placeholders like `INSTALL.md`, `USER_GUIDE.md`, and `CONTRIBUTING.md` with actual files or relevant sections in your repository. Additionally, include any specific dependencies or prerequisites required to run the software. 
