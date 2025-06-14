""" 
Automated image analysis pipeline (for macOS and Windows) for threshold-based quantification of 
pixels in images - for example fluorescence images. This script processes grayscale images by 
applying a binary threshold, measures fluorescent areas, saves the results and binary masks. It 
performs group-based statistical analysis and saves the data to Excel and creates a bar plot 
comparing group means and standard deviations. Pixel number, pixel size and image size must be the 
same for each image. The scaffold and code snippets of the script were generated with the assistance 
of AI and modified, extended and customized for error handling. As part of building this project, I 
created an annotated line by line walkthrough of the code to solidify my understanding.
"""

# Import standard and third-party libraries
import cv2
import numpy as np
import os
from glob import glob
import pandas as pd
from scipy.stats import shapiro, f_oneway, kruskal, levene, mannwhitneyu
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import statsmodels.stats.multitest as multitest
import matplotlib.pyplot as plt 
from itertools import combinations
 

# -------- Image processing --------

# Saves binary masks and returns DataFrame
def process_gfp_images(input_folder, output_folder, threshold=60, pixel_size_cm = 0.009025):

    # Create output folder if it does not exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Check if the folder was created 
    if not os.path.exists(output_folder):
        print(f"Failed to create output folder: {output_folder}")
    else:
        print(f"Output folder: {output_folder}")

    # Get all image paths
    image_paths = glob(os.path.join(input_folder, "*.jpg")) 
    print(f"Files found: {image_paths}")  
    if not image_paths:
        print(f"No images found in {input_folder}. Check the folder path.")
        return

    # Initialize a list to store results as dictionary
    results = [] 
    failed_images  = []

    # Loop through each image 
    for img_path in image_paths:
        print(f"Reading image: {img_path}") 
        
        image = cv2.imread(img_path, cv2.IMREAD_GRAYSCALE)

        if image is None:
            print(f"Failed to read {img_path}")
            failed_images.append(os.path.basename(img_path))
            continue                             

        # Apply threshold to create binary mask, pixels above threhold will have value of 255
        _, binary_mask = cv2.threshold(image, threshold, 255, cv2.THRESH_BINARY)
        
        # Save the binary mask image
        mask_filename = os.path.join(output_folder, os.path.basename(img_path).replace(".jpg", "_mask.jpg"))
        if cv2.imwrite(mask_filename, binary_mask):
            print(f"Saved mask: {mask_filename}") 

            cv2.imshow("Binary mask", binary_mask)
            cv2.waitKey(2000)
            cv2.destroyAllWindows() 

        else:
            print(f"Failed to save mask: {mask_filename}")  

        # Pixel area
        pixel_area_cm2 = pixel_size_cm ** 2 

        # Calculate fluorescent area as total area of white pixels in cm2
        white_pixel_count = np.sum(binary_mask == 255)
        white_pixel_area = white_pixel_count * pixel_area_cm2

        print(f"Processed {os.path.basename(img_path)} - Fluorescent Area: {white_pixel_count} \
            white pixels and {white_pixel_area:.2f} cm\u00B2.")


        # Append results (image name and fluorescent area) to the results list
        results.append({
            "Image Name": os.path.splitext(os.path.basename(img_path))[0], # Returns a tuple of root & extension 
            "Fluorescent Area cm\u00B2": white_pixel_area,
            "Fluorescent Area (pixel count)": white_pixel_count, 
        })

    if failed_images:
        df_failed = pd.DataFrame(failed_images, columns=["Failed image names"])
        failed_excel_path = os.path.join(output_folder, "failed_images.xlsx")
        df_failed.to_excel(failed_excel_path, index=False)
        print(f"Failed images saved to {failed_excel_path}") 
    

    # Convert results list to a pandas DataFrame
    df = pd.DataFrame(results)  
    print("Columns in DataFrame:", df.columns)
    return df



# -------- Statistical Analysis --------

def analyze_groups(df): 
    # Add group labels based on image name
    def get_group(name):
        if "_p" in name:
            return "Positive Control"
        elif "_n" in name:
            return "Negative Control"
        elif "_e" in name:
            return "Experiment"
        else:
            return "Unknown"

    df["Group"] = df["Image Name"].apply(get_group) 

    # Prepare stats collection
    group_stats = []
    normality_results = {}
    all_groups = df["Group"].unique()

    # Calculate stats, test normality and variance
    for group in all_groups:
        values = df[df["Group"] == group]["Fluorescent Area cm\u00B2"]
        mean = values.mean()
        std_dev = values.std()
        variance = values.var()
        
        # Shapiro-Wilk Test for Normality
        if len(values) >= 5:             # Min number should be adjusted according to requirements
            stat, p_value = shapiro(values) 
        else:                            # Fail fast and no misleading results; could also do graceful handling
            raise ValueError(f" Too few samples: {len(values)} for Shapiro-Wilk test in group '{group}."
                             f" At least 5 samples per group required.")

        group_stats.append({       
            'Group': group,
            'Mean': mean,
            'Standard Deviation': std_dev,
            'Variance': variance,
            'Shapiro-Wilk p-value': p_value
        })
        
        normality_results[group] = p_value

    # Test between-group differences
    test_result = {}
    grouped_values = [df[df['Group'] == group]['Fluorescent Area cm\u00B2'].values for group in all_groups]

    # Check if all groups are normally distributed (p > 0.05)
    normal = all(not np.isnan(p) and p > 0.05 for p in normality_results.values())

    # Levene's test for equal variances
    if len(values) >= 5:               # Min number should be adjusted according to requirements
        stat, levene_p = levene(*grouped_values)    
        equal_variance = levene_p > 0.05
    else:
        raise ValueError(f"Too few samples: {len(values)} for Levene's test in {group}. "
                        f"At least 5 samples per group required.")
    
    # ANOVA or Kruskal-Wallis test and Post-hoc analysis
    if normal and equal_variance:
        stat, p = f_oneway(*grouped_values)
        test_result = {'Test': 'ANOVA', 'Statistic': stat, 'p-value': p} 

        # Post-hoc Tukey HSD
        all_values = df[['Group', 'Fluorescent Area cm\u00B2']]
        tukey = pairwise_tukeyhsd(endog=all_values['Fluorescent Area cm\u00B2'], groups=all_values['Group'], alpha=0.05)
        
        # Extract results to DataFrame
        posthoc_df = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
    else:
        stat, p = kruskal(*grouped_values)
        test_result = {'Test': 'Kruskal-Wallis', 'Statistic': stat, 'p-value': p}

        # Pairwise Mann-Whitney U tests with Bonferroni correction
        pairs = []
        pvals = []
        groups = all_groups.tolist()
        for i in range(len(groups)):
            for j in range(i+1, len(groups)):
                group1 = groups[i]
                group2 = groups[j]
                vals1 = df[df['Group'] == group1]['Fluorescent Area cm\u00B2']
                vals2 = df[df['Group'] == group2]['Fluorescent Area cm\u00B2']
                stat_, p_ = mannwhitneyu(vals1, vals2, alternative='two-sided')
                pairs.append(f"{group1} vs {group2}")
                pvals.append(p_)

        # Multiple testing correction
        reject, pvals_corrected, _, _ = multitest.multipletests(pvals, alpha=0.05, method='bonferroni')
        posthoc_df = pd.DataFrame({
            'Comparison': pairs,
            'p-value raw': pvals,
            'p-value corrected': pvals_corrected,
            'Significant': reject
        })

    # Save results to Excel
    with pd.ExcelWriter('fluorescence_stats_output.xlsx') as writer:
        pd.DataFrame(group_stats).to_excel(writer, sheet_name='Group Stats', index=False)
        pd.DataFrame([test_result]).to_excel(writer, sheet_name='Group Comparison', index=False)
        df.to_excel(writer, sheet_name='Raw Data + Group', index=False)
        posthoc_df.to_excel(writer, sheet_name='Group Comparison Posthoc', index=False)
        print("Stats output file created.")
 

    # Convert to DataFrame
    stats_df = pd.DataFrame(group_stats)
    stats_df = stats_df.sort_values("Group") # Sort groups for consistent plotting orders



    # -------- Plotting the data ---------

    # Plotting with group specific colors
    group_colors = {
    'Positive Control': '#2ca02c',   
    'Negative Control': '#d62728',   
    'Experiment': '#1f77b4',         
    'Unknown': 'gray'
    }

    colors = stats_df["Group"].map(group_colors)  

    # Asymmetric error bars (lower = 0, upper = std dev)
    lower_errors = np.zeros_like(stats_df["Standard Deviation"].values)
    upper_errors = stats_df["Standard Deviation"].values
    asymmetric_errors = np.array([lower_errors, upper_errors])

    # Plot properties
    fig, ax = plt.subplots(figsize=(6, 6))
    bars = ax.bar(
        stats_df["Group"],
        stats_df["Mean"],
        yerr=asymmetric_errors,
        capsize=5,
        color= colors,
        edgecolor='black',
        ecolor='black',
        width=0.5
    )

    for i, bar in enumerate(bars):
        mean_val = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            -0.13 * max(stats_df["Mean"]),  # Negative offset
            f'{mean_val:.1f}',
            ha='center',
            va='top',  # 'top' so the text aligns above the baseline at y position
            fontsize=10,
        )

    # Add horizontal line at y = 0
    plt.axhline(0, color='gray', linewidth=1.5, linestyle='--')

    # Delete grid
    ax.grid(False) 

    # Labels and ticks properties
    ax.set_ylabel("Total Fluorescent Area/Image [cm\u00B2]", fontsize=14, fontweight='bold')
    ax.set_title("Distribution of Fluorescent Pixel Area by Group", fontsize=16)
    ax.tick_params(axis='x', bottom=False, top=False, labelsize=12)
    ax.tick_params(axis='y', length=4, labelsize=12)
    ax.tick_params(axis='y', length=4)

    for label in ax.get_xticklabels():
        label.set_fontweight('bold')

    # Optimize spacing of plot elements
    plt.tight_layout()

    # Process plot image
    plt.savefig("group_fluorescence_barplot.png", dpi=300)
    plt.show(block=False) # Code is not paused
    plt.pause(30)         # Seconds the plot window stays visible
    plt.close()


try: 
    if __name__ == "__main__":
        input_folder = "/Users/x/Documents/Coding/image_analysis/image_files"  
        output_folder = "./output_masks"

    # Process images
    df = process_gfp_images(input_folder, output_folder)
    
    # Confirm if file was created and contains rows, save and analyze results
    if not df.empty:
        df.to_excel("fluorescent_area_results.xlsx", index=False)   
        if not os.path.exists("fluorescent_area_results.xlsx"):
            print("Excel file not found after saving.")
            exit()
        print("Raw data saved to Excel file.")

        # Try reading the Excel file
        try:
            df_read = pd.read_excel("fluorescent_area_results.xlsx")
            print("Excel file successfully read.")

            if df_read.notna().sum().sum() == 0:
                print("Excel file has rows but all values are NaN.")
                exit()
            print(df_read.head())  # Optional: Check the first few rows of the DataFrame

        except Exception as e:
            print(f"Error reading the Excel file: {e}") 
            

        # Run statistical analysis
        analyze_groups(df)

    else:
        print("No image data processed. Nothing to save.") # If the DataFrame is empty
    
except Exception as e:
    print("\u274C An error occurred.")





