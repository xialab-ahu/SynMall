"""
@description: This script is used to assign a synScore value to all synonymous mutations in SynMall.
@author: chen ye
@date: 2025-10-03
"""


import os
import joblib
import numpy as np
import pandas as pd
import gzip
from tqdm import tqdm


def retrieve_variant(df, chrom_col, pos_col, ref_col, alt_col, prefix=''):
	"""
	@Desc:
		return a vid list: {chromosome}_{position}_{reference allele}/{{alternate allele}}
	@Args:
		df, dataframe
		chrom_col, name of chromosome column
		pos_col, name of position column
		ref_col, name of reference allele column
		alt_col, name of alternate allele column
	"""
	variant = [f'{prefix}{chro}_{pos}_{ref}/{alt}' for chro, pos, ref, alt in zip(df[chrom_col], df[pos_col], df[ref_col], df[alt_col])]
	return variant


def preprocess(df):
	"""
	@Desc:
		preprocess the dataframe as numeric features
	@Args:
		df, dataframe
	"""
	# 1 insert vid
	vid = retrieve_variant(df, '#CHROM', 'POS', 'REF', 'ALT')
	df.insert(0, 'Variant38', vid)
	# 2 get transcript info
	def parse_transcript_info(cell):
		if pd.isna(cell) or cell.strip() == "":
			return pd.Series({
				"d_esr": np.nan,
				"d_ess": np.nan,
				"d_ese": np.nan,
				"splice_site_acceptor": 0,
				"splice_site_donor": 0
			})
		entries = cell.split(";")
		best_ess, best_ese, best_esr = None, None, None
		splice_site_acceptor = 0
		splice_site_donor = 0
		for entry in entries:
			parts = entry.split("|")
			if len(parts) < 13:
				continue
			try:
				ess = float(parts[8])
				ese = float(parts[9])
				esr = float(parts[10])
				acceptor = parts[11].strip() == "True"
				donor = parts[12].strip() == "True"
				# Use the maximum value in all transcripts
				if best_ess is None or abs(ess) > abs(best_ess):
					best_ess = ess
				if best_ese is None or abs(ese) > abs(best_ese):
					best_ese = ese
				if best_esr is None or abs(esr) > abs(best_esr):
					best_esr = esr
				if acceptor:
					splice_site_acceptor = 1
				if donor:
					splice_site_donor = 1
			except Exception:
				continue
		return pd.Series({
			"d_esr": best_ess,
			"d_ess": best_ese,
			"d_ese": best_esr,
			"splice_site_acceptor": splice_site_acceptor,
			"splice_site_donor": splice_site_donor
		})
	if "transcript_wise_info" in df.columns:
		parsed = df["transcript_wise_info"].apply(parse_transcript_info)
		df = pd.concat([df, parsed], axis=1)
		df.drop(['#CHROM', 'POS', 'REF', 'ALT', 'transcript_wise_info'], axis=1, inplace=True)
	# 3 process non-numeric columns
	# super_enhancer: count items
	df['super_enhancer'] = df['super_enhancer'].apply(lambda x: 0 if pd.isna(x) else len(str(x).split(","))).astype(int)
	# genehancer: extract numeric score
	df['genehancer'] = pd.to_numeric(df['genehancer'].astype(str).str.extract(r"Score=([\d\.]+)")[0], errors='coerce')
	# cage_enhancer & cage_promoter: presence/absence
	df['cage_enhancer'] = (~df['cage_enhancer'].isna() & (df['cage_enhancer'].astype(str).str.strip() != '')).astype(int)
	df['cage_promoter'] = (~df['cage_promoter'].isna() & (df['cage_promoter'].astype(str).str.strip() != '')).astype(int)
	# ORegAnno_ID: count items separated by ;
	df['ORegAnno_ID'] = df['ORegAnno_ID'].apply(lambda x: 0 if pd.isna(x) else len(str(x).split(";"))).astype(int)
	df = df[sequence_header]
	return df


def rescale(df, imputer, scaler, empty_cols):
	"""
	@Desc:
		Impute and scale the original feature
	@Args:
		df, dataframe
		imputer, imputer of the training set
		scaler, scaler of the training set
		empty_cols, deleted columns in the training set because of all nan values.
	"""
	# split meta-info columns and feature columns
	df_info = df.iloc[:, :1]
	df_fea = df.iloc[:, 1:].copy()
	# replace nan value notation
	df_fea.replace(["na", "NA", "nan", ""], np.nan, inplace=True)
	# delete all empty cols (according to the previous training set)
	if empty_cols:
		df_fea.drop(columns=empty_cols, inplace=True, errors='ignore')
	# impute missing values based on training imputer
	df_fea_imputed = pd.DataFrame(imputer.transform(df_fea),
								  columns=df_fea.columns,
								  index=df_fea.index)
	# scale values based on training imputer
	df_fea_scaled = pd.DataFrame(scaler.transform(df_fea_imputed),
								 columns=df_fea_imputed.columns,
								 index=df_fea_imputed.index)
	# concatenate meta-info and preprocessed feature columns
	df_result = pd.concat([df_info, df_fea_scaled], axis=1)
	return df_result


def predict_ensemble_from_files(model_files, X_test):
	"""
	@Desc:
		Make predictions with an ensemble of LightGBM models loaded from files.
	@Args:
		model_files: list of str, paths to LightGBM model .pkl files
		X_test: pd.DataFrame, test features (columns must match training)
	@Returns:
		preds: averaged prediction probabilities
	"""
	models = [joblib.load(f) for f in model_files]
	preds = np.zeros((X_test.shape[0], len(models)))
	for i, model in enumerate(models):
		preds[:, i] = model.predict(X_test, num_iteration=getattr(model, "best_iteration", None))
	return preds.mean(axis=1)


if __name__ == '__main__':
	chr_list = ['chrY', 'chrX'] +  ['chr' + str(i) for i in range(1, 23)]
	input_dir = "../Result/tabix/"
	output_dir = "../Result/whole_genome/"
	os.makedirs(output_dir, exist_ok=True)

	# Read in columns
	original_header = pd.read_csv("../Data/MetaInfo/OriginHeader.txt", header=None)[0].tolist()
	use_header = pd.read_csv("../Data/MetaInfo/UseHeader.txt", header=None)[0].tolist()
	sequence_header = ['Variant38'] + pd.read_csv("../Data/MetaInfo/descriptors.txt", header=None)[0].tolist()

	# Load in models
	model_dir = "../Result/models/"
	model_files = [os.path.join(model_dir, f"lgb_fold{i}.pkl") for i in range(1, 6)]
	imputer = joblib.load("../Result/models/imputer.pkl")
	scaler = joblib.load("../Result/models/scaler.pkl")
	empty_cols = joblib.load("../Result/models/empty_cols.pkl")

	# Annotate data
	for chr_num in tqdm(chr_list):
		input_file = os.path.join(input_dir, f"{chr_num}.tsv.gz")
		output_file = os.path.join(output_dir, f"{chr_num}.synScore.txt")
		all_chunks = []
		# read in by chunks
		chunksize = 500000
		with gzip.open(input_file, "rt") as f:
			for chunk in pd.read_csv(f, sep='\t', header=None, comment='#', chunksize=chunksize, low_memory=False):
				# preprocess
				chunk.columns = original_header
				chunk = chunk[use_header]
				chunk_processed = preprocess(chunk)
				chunk_scaled = rescale(chunk_processed, imputer, scaler, empty_cols)
				# predict
				y_pred = predict_ensemble_from_files(model_files, chunk_scaled.iloc[:, 1:])
				chunk_processed['synScore'] = y_pred
				all_chunks.append(chunk_processed[['Variant38', 'synScore']])
		# Concatenate all results and save
		pd.concat(all_chunks, axis=0).to_csv(output_file, sep='\t', index=False)
		print(f"[INFO] Finished {chr_num}, saved to {output_file}")
