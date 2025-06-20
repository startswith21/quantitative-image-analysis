# Automated Image Analysis Pipeline for Signal Area Quantification

This project is a Python-based pipeline that analyzes grayscale images by applying a threshold to detect bright (white) pixels, quantifies white pixel areas and performs group-based statistical analysis with visual outputs. It works on macOS and Windows. Code logic can easily be adapted for other image types and measurements. Tested with Python 3.9.6.
Example use case: Measuring a signal in biological images to compare between control and experimental groups.


## Features
- Automatic processing of `.jpg` images
- Applies a binary threshold to extract bright pixels
- Calculates area in pixels and cm²
- Automatically classifies images into groups
- Performs statistical analysis with post-hoc tests
- Output: Binary images, Excel files and graph


## Requirements
Install dependencies using:
```bash
pip install numpy pandas opencv-python matplotlib scipy statsmodels
```


## Folder structure
```project-folder/
├── image_files/                     # Folder with input .jpg images 
├── output_masks/                    # Output folder for binary masks
├── fluorescent_area_results.xlsx
├── fluorescence_stats_output.xlsx
├── group_fluorescence_barplot.png
└── image_analysis_pipeline.py 
```

## Filename conventions
A group code must be included (_p, _n, _e):
- image01_p.jpg → Positive Control
- image02_n.jpg → Negative Control
- image03_e.jpg → Experimental Group


## Run the script
- Update the files`path:   
if __name__ == "__main__":
    input_folder = "/your/path/to/image_files"
    output_folder = "./output_masks".                  
- Run the script:    
python image_analysis_pipeline.py


## Pipeline process
### 1. Load images
Reads .jpg images from the image_files/ folder in grayscale

### 2. Apply threshold
Converts image to binary mask using a set threshold, which can be adjusted; *e.g.* 60:   
Pixel > 60 → white (255).     
Else → black (0)

### 3. Measure white pixels
Count and calculate total area; pixel size must be adjusted and all images must have the same pixel size

### 4. Save binary masks
Saves binary versions of each image in output_masks/ folder

### 5. Classify images by group based on filename
_p → Positive Control.  
_n → Negative Control.  
_e → Experiment.  

### 6. Statistical hypothesis testing
- Shapiro-Wilk test for normality
- Levene’s test for equal variance
- ANOVA or Kruskal-Wallis for testing null hypothesis
- Tukey HSD or Mann-Whitney U for post-hoc comparison
- Write results to Excel file.    
→ Code is currently set to a pragmatic lower limit of 5 images per group. Recommended to collect more 
  images, *i.e.* 8 - 10.

### 7. Plot the data
Bar plot comparing group means and standard deviation


## Credits
- Built with the help of OpenAI-generated scaffolds
- Customized, extended and debugged manually
- Created to better understand statistics, plotting and use of Excel with Python