
# Load libraries
import os
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from scipy.stats import zscore

# Set file path
data_dir   = r"C:\Users\Arthur\OneDrive\Documents\Python Scripts\LDT_Generator\Data"
output_dir = r"C:\Users\Arthur\OneDrive\Documents\Python Scripts\LDT_Generator\Output"

# Load data frames
elp_words     = pd.read_csv(os.path.join(data_dir, "ELP_words.csv"))
flp_words     = pd.read_csv(os.path.join(data_dir, "FLP_words.csv"))
# fr_aoa_data   = pd.read_csv(os.path.join(data_dir, "alarioFerrand1999_aoa.tsv"), sep="\t")
fr_concr_data = pd.read_csv(os.path.join(data_dir, "boninEtAl2003_concr_imag.csv"))
elp_nonwords  = pd.read_csv(os.path.join(data_dir, "ELP_nonwords.csv"))
flp_nonwords  = pd.read_csv(os.path.join(data_dir, "FLP_pseudowords.csv"))

# Reduce the ELP words df to nouns only
elp_words = elp_words[elp_words["POS"] == "NN"]

# Keep only necessary columns
elp_words = elp_words[["Word","Length","NSyll","LgSUBTLWF",
                       "I_Mean_RT","Age_Of_Acquisition","Concreteness_Rating"]]
flp_words = flp_words[["item","nletters","nsyllables","lcfreqmovies","rt"]]
# fr_aoa_data = fr_aoa_data[["Word","aa_m"]]
fr_concr_data = fr_concr_data[["Word","Concr.M"]]
elp_nonwords = elp_nonwords[["Word","Length","NWI_Mean_RT"]]
flp_nonwords = flp_nonwords[["item","rt"]]

# Rename cols
elp_words = elp_words.rename(columns={
    "Word":"word",
    "Length":"nletters",
    "NSyll":"nsyllables",
    "LgSUBTLWF":"log_freq",
    "I_Mean_RT":"rt_mean",
    # "Age_Of_Acquisition": "aoa",
    "Concreteness_Rating": "concr"
})
flp_words = flp_words.rename(columns={
    "item":"word",
    "lcfreqmovies":"log_freq",
    "rt":"rt_mean"
})
# fr_aoa_data = fr_aoa_data.rename(columns={
#     "Word":"word",
#     "aa_m":"aoa"
# })
fr_concr_data = fr_concr_data.rename(columns={
    "Word":"word",
    "Concr.M":"concr"
})
elp_nonwords = elp_nonwords.rename(columns={
    "Word":"word",
    "Length":"nletters",
    "NWI_Mean_RT":"rt_mean"
})
flp_nonwords = flp_nonwords.rename(columns={
    "item":"word",
    "rt":"rt_mean"
})

# Compute length for flp_nonwords, since it is missing
flp_nonwords["nletters"] = flp_nonwords["word"].str.len()
flp_nonwords = flp_nonwords[["word","nletters","rt_mean"]]

# Function to normalize words to prevent merge failures
def normalize_word(word):

    # Convert word to lowercase and strip extra white space
    word = word.lower().strip()

    # Was going to also standardize accents but wasn't necessary

    return word

# Apply the above function to the 3 dfs of French words
flp_words["word"]     = flp_words["word"].apply(normalize_word)
# fr_aoa_data["word"]   = fr_aoa_data["word"].apply(normalize_word)
fr_concr_data["word"] = fr_concr_data["word"].apply(normalize_word)

# Perform inner join on the 3 dfs of French words
# fr_merged_df = pd.merge(flp_words,    fr_aoa_data,   on="word", how="inner")
# fr_merged_df = pd.merge(fr_merged_df, fr_concr_data, on="word", how="inner")
fr_merged_df = pd.merge(flp_words, fr_concr_data, on="word", how="inner")

# Ensure that all cols except "word" are numeric
cols_to_convert = elp_words.columns.drop("word")
elp_words[cols_to_convert] = elp_words[cols_to_convert].apply(
    pd.to_numeric, errors="coerce")
cols_to_convert = fr_merged_df.columns.drop("word")
fr_merged_df[cols_to_convert] = fr_merged_df[cols_to_convert].apply(
    pd.to_numeric, errors="coerce")
cols_to_convert = elp_nonwords.columns.drop("word")
elp_nonwords[cols_to_convert] = elp_nonwords[cols_to_convert].apply(
    pd.to_numeric, errors="coerce")
cols_to_convert = flp_nonwords.columns.drop("word")
flp_nonwords[cols_to_convert] = flp_nonwords[cols_to_convert].apply(
    pd.to_numeric, errors="coerce")

# Convert infinite values to NaNs
elp_words["log_freq"] = (
    elp_words["log_freq"]
    .replace([np.inf, -np.inf], np.nan))
fr_merged_df["log_freq"] = (
    fr_merged_df["log_freq"]
    .replace([np.inf, -np.inf], np.nan))

# Drop rows with NaN values
elp_words    = elp_words.dropna()
fr_merged_df = fr_merged_df.dropna()
elp_nonwords = elp_nonwords.dropna()
flp_nonwords = flp_nonwords.dropna()

# Compare the attributes of the two finished dataframes
print(elp_words["rt_mean"].mean())
print(fr_merged_df["rt_mean"].mean())
print(elp_words["log_freq"].mean())
print(fr_merged_df["log_freq"].mean())

# Log frequency is evidently not comparable b/w the two datasets, so switching to z-scores
elp_words["log_freq_z"]    = zscore(elp_words["log_freq"])
fr_merged_df["log_freq_z"] = zscore(fr_merged_df["log_freq"])

# Create thresholds for high- and low-frequency words
THRESH_EN_HIGH = elp_words["log_freq_z"].quantile(.70)
THRESH_EN_LOW  = elp_words["log_freq_z"].quantile(.30)
THRESH_FR_HIGH = fr_merged_df["log_freq_z"].quantile(.70)
THRESH_FR_LOW  = fr_merged_df["log_freq_z"].quantile(.30)

# Create high- and low-frequency dfs
en_hf_pool = elp_words[elp_words["log_freq_z"] >= THRESH_EN_HIGH]
en_lf_pool = elp_words[elp_words["log_freq_z"] <= THRESH_EN_LOW]
fr_hf_pool = fr_merged_df[fr_merged_df["log_freq_z"] >= THRESH_FR_HIGH]
fr_lf_pool = fr_merged_df[fr_merged_df["log_freq_z"] <= THRESH_FR_LOW]

# Function to created matched lists of stimuli based on Mahalanobis distance
def mahalanobis_match(df1, df2, match_vars, n):

    # Extract the variables for matching and convert them to a NumPy matrix
    x1 = df1[match_vars].to_numpy()
    x2 = df2[match_vars].to_numpy()

    # Combine the two matrices
    x_all = np.vstack([x1, x2])
    # Compute the covariance matrix
    cov = np.cov(x_all, rowvar=False)
    # Compute the inverse covariance matrix
    inv_cov = np.linalg.pinv(cov)

    # Calculate the Mahalanobis distance between every possible word pair
    dist_matrix = cdist(
        x1,
        x2,
        metric="mahalanobis",
        VI=inv_cov
    )

    # Obtain optimal one-to-one matches
    row_ind, col_ind = linear_sum_assignment(dist_matrix)
    pairs = pd.DataFrame({
        "df1_row": row_ind,
        "df2_row": col_ind,
        "distance": dist_matrix[row_ind, col_ind]
    })
    pairs = pairs.sort_values("distance").head(n)

    # Extract the matched words from the original dataframes
    matched_1 = df1.iloc[pairs["df1_row"]].copy()
    matched_2 = df2.iloc[pairs["df2_row"]].copy()
    # Assign match_id to show which words go together
    matched_1["match_id"] = range(1, len(matched_1) + 1)
    matched_2["match_id"] = range(1, len(matched_2) + 1)

    return matched_1, matched_2, pairs

# Extract matching high-freq English words and high-freq French words
en_hf_matched, fr_hf_matched, hf_pairs = mahalanobis_match(
    en_hf_pool,
    fr_hf_pool,
    match_vars=["rt_mean","log_freq_z","nletters","nsyllables","concr"],
    n=32
)

# Extract matching low-freq English words and low-freq French words
en_lf_matched, fr_lf_matched, lf_pairs = mahalanobis_match(
    en_lf_pool,
    fr_lf_pool,
    match_vars=["rt_mean","log_freq_z","nletters","nsyllables","concr"],
    n=32
)

# Create one df with all words, with cols for language and condition
final_words = pd.concat([
    en_hf_matched.assign(language="English", condition="HF"),
    fr_hf_matched.assign(language="French",  condition="HF"),
    en_lf_matched.assign(language="English", condition="LF"),
    fr_lf_matched.assign(language="French",  condition="LF")
])

# Print the means and SDs by condition for the matching variables
pd.set_option("display.max_columns", None)
print(final_words.groupby(["language", "condition"])[
    ["rt_mean","log_freq_z","nletters","nsyllables","concr"]
].agg(["mean", "std"]))
pd.reset_option("display.max_columns")

# Select a subset of English nonwords at random
# Because the function linear_sum_assignment() was prohibitively slow with both full datasets of nonwords
# random_state is setting a seed to ensure reproducibility
elp_nw_sample = elp_nonwords.sample(n=500, replace=False, random_state=42)

# Extract matching "English" nonwords and "French" nonwords
en_nw_matched, fr_nw_matched, nw_pairs = mahalanobis_match(
    elp_nw_sample,
    flp_nonwords,
    match_vars=["rt_mean","nletters"],
    n=32
)

# Create one df with all nonwords, with col for language
final_nonwords = pd.concat([
    elp_nonwords.assign(language="English"),
    flp_nonwords.assign(language="French")
])

# Print the means and SDs by condition for the matching variables
pd.set_option("display.max_columns", None)
print(final_nonwords.groupby(["language"])[
    ["rt_mean","nletters"]
].agg(["mean", "std"]))
pd.reset_option("display.max_columns")

# Export the final word lists as .csv files
# Note that there are 32 per condition - this is in case any words need to be dropped (otherwise the last 2 will be dropped)
en_hf_matched.to_csv(os.path.join(output_dir, "english_high_frequency_words.csv"),
                     index=False, encoding="utf-8-sig")
en_lf_matched.to_csv(os.path.join(output_dir, "english_low_frequency_words.csv"),
                     index=False, encoding="utf-8-sig")
fr_hf_matched.to_csv(os.path.join(output_dir, "french_high_frequency_words.csv"),
                     index=False, encoding="utf-8-sig")
fr_lf_matched.to_csv(os.path.join(output_dir, "french_low_frequency_words.csv"),
                     index=False, encoding="utf-8-sig")
en_nw_matched.to_csv(os.path.join(output_dir, "english_nonwords.csv"),
                     index=False, encoding="utf-8-sig")
fr_nw_matched.to_csv(os.path.join(output_dir, "french_nonwords.csv"),
                     index=False, encoding="utf-8-sig")


