""" 
Automated image analysis pipeline (for macOS and Windows) for threshold-based quantification of 
fluorescence signals in microscopy images, processing datasets of single-channel fluorescence images. 
This script processes grayscale fluorescence images by applying a binary threshold,
measures fluorescent areas, saves the results and binary masks. It performs group-based
statistical analysis (normality, ANOVA/Kruskal-Wallis) and saves the data to Excel 
and creates a bar plot comparing group means and standard deviations. The scaffold and code snippets 
of the script were generated with the assistance of AI and modified, extended and customized.
As part of building this project, I created an annotated line by line walkthrough of the 
code to solidify my understanding.
"""

import cv2
import numpy as np
import os
from glob import glob
import pandas as pd
from scipy.stats import shapiro, f_oneway, kruskal
import matplotlib.pyplot as plt 


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
            "Image Name": os.path.splitext(os.path.basename(img_path))[0], # returns a tuple of root & extension 
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

    # Calculate stats and test normality
    for group in all_groups:
        values = df[df["Group"] == group]["Fluorescent Area (pixel count)"]
        mean = values.mean()
        std_dev = values.std()
        variance = values.var()
        
        # Shapiro-Wilk Test for Normality
        if len(values) >= 3:
            stat, p_value = shapiro(values)
        else:
            stat, p_value = np.nan, np.nan  # too few samples
            print("Too few samples.")
            exit()

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
    grouped_values = [df[df['Group'] == group]['Fluorescent Area (pixel count)'].values for group in all_groups]

    # Check if all groups are normal (p > 0.05)
    normal = all(p is not np.nan and p > 0.05 for p in normality_results.values())

    if normal:
        stat, p = f_oneway(*grouped_values)
        test_result = {'Test': 'ANOVA', 'Statistic': stat, 'p-value': p} 
    else:
        stat, p = kruskal(*grouped_values)
        test_result = {'Test': 'Kruskal-Wallis', 'Statistic': stat, 'p-value': p}

    # Save results to Excel
    with pd.ExcelWriter('fluorescence_stats_output.xlsx') as writer:
        pd.DataFrame(group_stats).to_excel(writer, sheet_name='Group Stats', index=False)
        pd.DataFrame([test_result]).to_excel(writer, sheet_name='Group Comparison', index=False)
        df.to_excel(writer, sheet_name='Raw Data + Group', index=False)
        print("Stat excel created")


    # Convert to DataFrame
    stats_df = pd.DataFrame(group_stats)
    stats_df = stats_df.sort_values("Group") # Sort groups for consistent plotting orders

    # Plot properties
    plt.figure(figsize=(8, 6))
    bars = plt.bar(
        stats_df["Group"],
        stats_df["Mean"],
        yerr=stats_df["Standard Deviation"],
        capsize=5,
        color='gray',
        edgecolor='black',
        width=0.4  
    ) 
    # Show mean vaue in graph
    for bar in bars:
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            height,
            f'{height:.1f}',
            ha='center',
            va='bottom',
            fontsize=10
    )
    # Add horizontal line at y = 0
    plt.axhline(0, color='gray', linewidth=1.5, linestyle='--')

    # Labels with font sizes
    plt.xlabel("Group", fontsize=14)
    plt.ylabel("Mean Fluorescent Area (pixel count)", fontsize=14)
    plt.title("Mean ± Standard Deviation of Fluorescent Area by Group", fontsize=16)

    # Bigger tick labels
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    #plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig("group_fluorescence_barplot.png", dpi=300)

    plt.show(block=False) # Code is not paused
    plt.pause(30)  # Seconds the plot window stays visible
    plt.close()


if __name__ == "__main__":
    input_folder = "/Users/x/Documents/Coding/image_analysis/image_files"  
    output_folder = "./output_masks"
    
    # Process images
    df = process_gfp_images(input_folder, output_folder) # Threshold can be overwritten by adding: ,70
    
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
            exit()
        
        # Run statistical analysis
        analyze_groups(df)

    else:
        # If the DataFrame is empty
        print("No image data processed. Nothing to save.")




