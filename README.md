## Paper code

Code for all papers is kept within the namespace for that paper as `+Author_Year`

To generate figures: 
- Download data
  - Download from relevant OSF project as `Author Year` at this address: https://osf.io/user/pkh5w
  - Data in this folder was packaged using the code namespace `+Author_Year.package` (this does not work without access to raw data)
- Run script which generates all figures in each paper: `+Author_Year.figures.generate_figures`
  - Data paths are usually sest at the top of this script, and will need to be changed to point to locally downloaded data
  - Code may have external dependencies, which are listed in the comments at the top of the figure script
