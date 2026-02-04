
#!/usr/bin/env python3
"""
Classifier.py  (FULL Python3 update + "pipeline-agnostic" table mode)

This file is a Python3-compatible update of the original Levitan/Berman-lab

  write_auc_curve, write_ortholog_excel, main, COLS_CONFIG, ORTH_COLS_CONFIG).
- All Python2-only syntax is fixed (print, xrange, iteritems, itervalues, tuple
  unpacking in def, lambda arg unpacking, cPickle, sklearn.cross_validation).
- A new **table mode** is added so you can:
    * Train on S. cerevisiae (dependencies / Organisms) and
    * Classify one or more **C. glabrata BG2 feature tables** (e.g. *.feature_table.RDF_1.csv)
    * Output per-library tables.xlsx AND optionally a combined table across libraries
    * Optionally join BG2↔Scer orthology (your C_glabrata_BG2_S_cerevisiae_orthologs.txt)

USAGE 
------------------------------------------------------
python3 Classifier.py \
  --mode scer-train-gla-classify \
  --gla-feature-glob "/path/to/04.summary/*/*.feature_table.RDF_1.csv" \
  --out-dir "/path/to/06.classifier" \
  --threads 32 \
  --target-fpr 0.10 \
  --combine \
  --orthology-file "/path/to/C_glabrata_BG2_S_cerevisiae_orthologs.txt"
  --overwrite \

USAGE in legacy paper mode; runs the original main() logic
--------------------------------------------------------
python3 Classifier.py --mode legacy-paper

NOTES
-----
1) Table mode expects glabrata feature tables with columns at least:
   standard_name, coding_length, hits, reads, neighborhood_index,
   freedom_index, upstream_hits_100
   (If coding_length is present, we map it to "length" used by feature groups.)

2) We keep the original "G4" group as default for table mode:
     ("neighborhood_index","length","hits","reads","freedom_index","upstream_hits_100")

3) Threads are applied via RandomForestClassifier(n_jobs=threads).
"""

import Shared
import GenomicFeatures
import SummaryTable
import Organisms

import glob
import os
import csv
import shutil
import re
from itertools import product, chain, repeat
from collections import OrderedDict
from pprint import pprint

# ---- Python2 -> Python3: cPickle ----
try:
    import cPickle as pickle  # type: ignore
except Exception:
    import pickle  # noqa: F401

# ---- sklearn: old vs new StratifiedKFold ----
try:
    from sklearn.model_selection import StratifiedKFold
    _SKF_STYLE = "new"
except Exception:
    from sklearn.cross_validation import StratifiedKFold  # type: ignore
    _SKF_STYLE = "old"

from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib_venn import venn2, venn3  # noqa: E402


# ----------------------------
# Python2 compatibility shims
# ----------------------------
def _iteritems(d):
    return getattr(d, "iteritems", d.items)()


def _itervalues(d):
    return getattr(d, "itervalues", d.values)()


def _xrange(n):
    return range(n)


# ----------------------------
# Core classifier helpers
# ----------------------------
def test_classifier(classifier, features, annotations):
    """Test a classifier using given features and annotations. Uses a 5-fold
    cross-validation method to score all of the test data points.

    Parameters
    ----------
    classifier : sklearn Classifier
        The classifier to test.
    features : list of sequences
        A list of features per gene.
    annotations : list of 1 or 0
        A list of annotations for the genes.

    Returns
    -------
    tuple of tuples
        The ROC curve parameters, as returned by `sklearn.metrics.roc_curve`,
        and the scores, as the second item of the tuple.
    """
    scores = np.empty(len(annotations), dtype=float)

    if _SKF_STYLE == "new":
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
        splits = skf.split(np.zeros(len(annotations)), annotations)
    else:
        skf = StratifiedKFold(y=annotations, n_folds=5, shuffle=True, random_state=0)  # type: ignore
        splits = skf  # old style yields (train_ix, test_ix)

    for roc_train, roc_test in splits:
        roc_train_set = [features[i] for i in roc_train]
        roc_test_set = [features[i] for i in roc_test]
        roc_train_essentiality = [annotations[i] for i in roc_train]

        classifier.fit(roc_train_set, roc_train_essentiality)
        scores[roc_test] = classifier.predict_proba(roc_test_set)[:, 1]

    return roc_curve(annotations, scores), scores


def get_training_groups():
    """Get all interesting permutations of features to train on.

    Note
    ----
    Use with caution, as this can be quite large.

    Returns
    -------
    generator of tuples
        Each tuple contains a permutation of features to use for training.
    """
    subgroups = (
        ("neighborhood_index", "insertion_index"),
        ("length", "hits", ("length", "hits"), ""),
        ("max_free_region", "freedom_index", "adj_max_free_region", "adj_freedom_index", ""),
        ("upstream_hits_50", "upstream_hits_100", ("upstream_hits_50", "upstream_hits_100"), ""),
        ("std_dev", "index_of_dispersion", ""),
    )

    def _flatten(features):
        # One-level flattening (sufficient for our structure).
        for f in features:
            if isinstance(f, str):
                yield f
            else:
                for sub_f in f:
                    yield sub_f

    for group in product(*subgroups):
        yield tuple(f for f in _flatten(group) if f)


def get_youden_statistic(fpr, tpr, thresholds):
    """The Youden index (or J statistic) is the threshold which maximizes
    the (TPR-FPR) value."""
    return max(zip(thresholds, tpr, fpr), key=lambda x: x[1] - x[2])


def draw_venn(title, output_file_path, data, labels):
    """A convenience function for generating Venn diagrams and saving them to
    a destination on the drive. The `data` and `labels` can be of size 2 or 3.

    Parameters
    ----------
    title : str
        The title of the plot.
    output_file_path : str
        The path in which the plot will be saved.
    data : list of set
        The data to plot - should be represented as sets.
    labels : list of str
        The labels of the data to plot.
    """
    plt.figure(figsize=(8, 8))

    venn = venn2 if len(data) == 2 else venn3
    v = venn(data, labels)

    for l in v.set_labels:
        l.set_style("italic")

    plt.title(title)
    plt.savefig(output_file_path, transparent=True, dpi=300)
    plt.close()


def write_group_legends(feature_groups, output_file):
    with open(output_file, "w", newline="") as out_file:
        writer = csv.writer(out_file, delimiter=",")
        writer.writerow(["Group name", "Features"])
        # feature_groups may be dict or list of tuples in some call paths
        if hasattr(feature_groups, "items"):
            it = feature_groups.items()
        else:
            it = feature_groups
        for group_name, columns_for_classification in it:
            writer.writerow([group_name, ", ".join(columns_for_classification)])


def only_orfs(records):
    return (r for r in records if "ORF" in r["feature"].type)


def no_dubious_features(records):
    return (r for r in records if getattr(r["feature"], "feature_qualifier", "") != "Dubious")


def no_depleted_hit_features(records):
    # If there were no insertions in the feature AND its neighborhood,
    # we can't trust this feature (maybe it was deleted):
    return (r for r in records if r["neighborhood_hits"] > 0 or r["hits"] > 0)


def combine_pre_post(pre_records, post_records):
    """Combine pre- and post- records into a single record list, by appending
    prefixes to the record values."""
    pre_records_sorted = sorted(pre_records, key=lambda r: r["feature"].standard_name)
    post_records_sorted = sorted(post_records, key=lambda r: r["feature"].standard_name)
    assert len(pre_records_sorted) == len(post_records_sorted)

    result = []
    for pre_record, post_record in zip(pre_records_sorted, post_records_sorted):
        assert pre_record["feature"] == post_record["feature"]
        new_record = {}
        for key, value in pre_record.items():
            new_record[key] = value
            new_record["pre_" + key] = value
        for key, value in post_record.items():
            new_record["post_" + key] = value
        result.append(new_record)

    return result


def explode_col_config(col_config):
    """"Explodes" the column configuration to display pre- and post- data, by
    appending prefixes to the column names."""
    result = []

    for config in col_config:
        if config.get("local", False):
            for prefix in ("pre", "post"):
                new_config = dict(config)
                new_config["field_name"] = "%s_%s" % (prefix, new_config["field_name"])
                new_config["csv_name"] = "%s - %s" % (prefix.title(), new_config["csv_name"])
                result.append(new_config)
        else:
            result.append(config)

    return result


def run_pipeline(cls_factory, cls_feature_groups, data, train_col_config, class_col_config, output_folder, target_fpr):
    """Classify the given datasets and store the output into the given folder.

    Parameters
    ----------
    cls_factory : dict of str to factory function
    cls_feature_groups : dict of str to tuple
    data : dict of str to dataset descriptor
    train_col_config : list[dict]
    class_col_config : list[dict]
    output_folder : str
    target_fpr : float
    """
    # Insert the classifier data into the column configuration.
    # We assume that the last two columns are type and description, so we want them to be last.
    train_col_config = list(train_col_config)
    class_col_config = list(class_col_config)

    train_col_config[-2:-2] = (
        [{"field_name": "ground_truth", "csv_name": "Ess. ground truth"}]
        + sum(
            [
                [
                    {"field_name": "train-%s-%s" % key, "csv_name": "%s - %s" % key, "format": "%.3f", "local": True},
                    {
                        "field_name": "train-%s-%s-verdict" % key,
                        "csv_name": "%s - %s - ess. for FPR %.3f" % (key + (target_fpr,)),
                        "local": True,
                    },
                ]
                for key in product(cls_factory.keys(), cls_feature_groups.keys())
            ],
            [],
        )
    )

    class_col_config[-2:-2] = (
        [{"field_name": "ground_truth", "csv_name": "Ess. ground truth"}]
        + sum(
            [
                [
                    {"field_name": "%s-%s" % key, "csv_name": "%s - %s" % key, "format": "%.3f", "local": True},
                    {
                        "field_name": "%s-%s-verdict" % key,
                        "csv_name": "%s - %s - ess. for FPR %.3f" % (key + (target_fpr,)),
                        "local": True,
                    },
                ]
                for key in product(cls_factory.keys(), cls_feature_groups.keys())
            ],
            [],
        )
    )

    threshold_ixs = {}
    for data_key in data.keys():
        data_folder = os.path.join(output_folder, "classification - %s" % data_key)
        Shared.make_dir(data_folder)

        train_ess_records, train_non_ess_records, class_datasets, benchmarks = (data[data_key] + ({},))[:4]
        train_all_records = train_ess_records + train_non_ess_records
        train_annotations = [1] * len(train_ess_records) + [0] * len(train_non_ess_records)

        for grp_name, features in cls_feature_groups.items():
            train_all_features = [[r[f] for f in features] for r in train_all_records]

            for cls_name, cls_creator in cls_factory.items():
                record_label = "%s-%s" % (cls_name, grp_name)
                train_record_label = "train-" + record_label
                verdict_record_label = record_label + "-verdict"
                train_verdict_record_label = "train-" + record_label + "-verdict"

                # Test classifier on training data
                ((fprs, tprs, thresholds), scores) = test_classifier(
                    cls_creator(),
                    train_all_features,
                    train_annotations,
                )

                # Write out AUC curve
                write_auc_curve(
                    (fprs, tprs, thresholds),
                    os.path.join(data_folder, "train_AUC.%s.png" % record_label),
                )

                training_threshold_ix = np.absolute(fprs - target_fpr).argmin()
                threshold_ixs[(cls_name, grp_name)] = training_threshold_ix
                fpr_threshold = thresholds[training_threshold_ix]

                for r, s in zip(train_all_records, scores):
                    r[train_record_label] = s
                    r[train_verdict_record_label] = "Yes" if s >= fpr_threshold else "No"

                # Train classifier
                classifier = cls_creator()
                classifier.fit(train_all_features, train_annotations)

                # Feature importances?
                if hasattr(classifier, "feature_importances_"):
                    print(data_key, cls_name, grp_name)
                    pprint(dict(zip(features, classifier.feature_importances_)))
                    print("\n")
                elif hasattr(classifier, "coef_"):
                    print(data_key, cls_name, grp_name)
                    pprint(dict(zip(features, classifier.coef_[0])))
                    print("\n")

                # Classify records for classification, store in records
                for classification_records in chain(*class_datasets.values()):
                    classification_features = [[r[f] for f in features] for r in classification_records]
                    predictions = classifier.predict_proba(classification_features)
                    for record, prediction in zip(classification_records, predictions[:, 1]):
                        record[record_label] = prediction
                        record[verdict_record_label] = "Yes" if prediction >= fpr_threshold else "No"

                for bench_label, (bench_ess, bench_non_ess) in benchmarks.items():
                    bench_all_records = bench_ess + bench_non_ess
                    bench_all_features = [[r[f] for f in features] for r in bench_all_records]
                    predictions = classifier.predict_proba(bench_all_features)
                    for record, prediction in zip(bench_all_records, predictions[:, 1]):
                        record[record_label] = prediction
                        record[verdict_record_label] = "Yes" if prediction >= fpr_threshold else "No"

                    roc_data = roc_curve(
                        [1] * len(bench_ess) + [0] * len(bench_non_ess),
                        predictions[:, 1],
                    )

                    write_auc_curve(
                        roc_data,
                        os.path.join(data_folder, "bench_AUC.%s.%s.png" % (bench_label, record_label)),
                    )

        with pd.ExcelWriter(os.path.join(data_folder, "tables.xlsx")) as excel_writer:
            # In Excel, the order of the sheets is:
            # 1) Prediction table - combined.
            # 2) Prediction table - pre.
            # 3) Prediction table - post.
            # 4) Training table
            # 5) FPR tables - classifiers X feature groups

            combined_col_config = explode_col_config(class_col_config)
            combined_col_config[-2:-2] = [
                {
                    "field_name": "%s-%s-verdict" % key,
                    "csv_name": "%s - %s - final ess. verdict" % key,
                }
                for key in product(cls_factory.keys(), cls_feature_groups.keys())
            ]

            # TODO: can't be generalized at all! (kept as-is)
            for ds_label, datasets in class_datasets.items():
                if len(datasets) == 2:
                    combined_records = combine_pre_post(datasets[0], datasets[1])
                    for (cls_name, grp_name) in product(cls_factory.keys(), cls_feature_groups.keys()):
                        pre_verdict_key = "pre_%s-%s-verdict" % (cls_name, grp_name)
                        post_verdict_key = "post_%s-%s-verdict" % (cls_name, grp_name)
                        combined_verdict_key = "%s-%s-verdict" % (cls_name, grp_name)
                        for r in combined_records:
                            combined_verdict = {
                                ("Yes", "Yes"): "Essential",
                                ("No", "No"): "Not essential",
                                ("Yes", "No"): "Weirdo",
                                ("No", "Yes"): "Sick",
                            }[(r[pre_verdict_key], r[post_verdict_key])]
                            r[combined_verdict_key] = combined_verdict
                    combined_df = SummaryTable.write_data_to_data_frame(combined_records, combined_col_config)
                    combined_df.to_excel(excel_writer, sheet_name="%s - combined class." % ds_label, index=False)
                else:
                    assert len(datasets) == 1
                    pred_df = SummaryTable.write_data_to_data_frame(datasets[0], class_col_config)
                    pred_df.to_excel(excel_writer, sheet_name="%s - class" % ds_label, index=False)

            train_df = SummaryTable.write_data_to_data_frame(train_ess_records + train_non_ess_records, train_col_config)
            train_df.to_excel(excel_writer, sheet_name="Training", index=False)

            for cls_name, grp_name in product(cls_factory.keys(), cls_feature_groups.keys()):
                label = "train-%s-%s" % (cls_name, grp_name)
                fprs, tprs, thresholds = roc_curve(
                    train_annotations,
                    [r[label] for r in train_all_records],
                )
                was_selected = [""] * len(fprs)
                was_selected[threshold_ixs[(cls_name, grp_name)]] = "Selected"
                fpr_df = pd.DataFrame(
                    {"FPR": fprs, "TPR": tprs, "Threshold": thresholds, "Was selected": was_selected},
                    columns=["FPR", "TPR", "Threshold", "Was selected"],
                )
                fpr_df.to_excel(excel_writer, sheet_name="FPRs - %s - %s" % (cls_name, grp_name), index=False)

            # Print benchmark tables
            for bench_label, (bench_ess, bench_non_ess) in benchmarks.items():
                bench_df = SummaryTable.write_data_to_data_frame(bench_ess + bench_non_ess, class_col_config)
                bench_df.to_excel(excel_writer, sheet_name="Benchmark %s - class." % bench_label, index=False)

                bench_all_records = bench_ess + bench_non_ess
                bench_annotations = [1] * len(bench_ess) + [0] * len(bench_non_ess)

                for cls_name, grp_name in product(cls_factory.keys(), cls_feature_groups.keys()):
                    label = "%s-%s" % (cls_name, grp_name)
                    fprs, tprs, thresholds = roc_curve(
                        bench_annotations,
                        [r[label] for r in bench_all_records],
                    )
                    fpr_df = pd.DataFrame(
                        {"FPR": fprs, "TPR": tprs, "Threshold": thresholds},
                        columns=["FPR", "TPR", "Threshold"],
                    )
                    fpr_df.to_excel(
                        excel_writer,
                        sheet_name="%s - FPRs-%s-%s" % (bench_label, cls_name, grp_name),
                        index=False,
                    )


def write_auc_curve(roc_tuple, output_file, title=None):
    (fpr, tpr, _threshold) = roc_tuple
    classifier_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, label="AUC = %.3f" % classifier_auc)
    plt.legend(loc="lower right")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")

    if title is not None:
        plt.title(title)

    plt.savefig(output_file, transparent=True, dpi=300)
    plt.close()


def write_ortholog_excel(orth_df, calb_fprs_df, scer_fprs_df, spom_fprs_df, output_file):
    orth_sheet_name = "Orthologs"
    ctrl_sheet_name = "Controls"
    calb_fprs_sheet_name = "Calb FPRs"
    scer_fprs_sheet_name = "Scer FPRs"
    spom_fprs_sheet_name = "Spom FPRs"

    with pd.ExcelWriter(output_file) as writer:
        workbook = writer.book
        for prefix in ("Ca", "Sc", "Sp"):
            insertion_index = list(orth_df.columns.values).index("%s RF G4" % prefix) + 1
            orth_df.insert(loc=insertion_index, column="%s verdict" % prefix, value=0)
        orth_df.to_excel(writer, sheet_name=orth_sheet_name, index=False)
        ctrl_sheet = workbook.add_worksheet(ctrl_sheet_name)
        calb_fprs_df.to_excel(writer, sheet_name=calb_fprs_sheet_name, index=False)
        scer_fprs_df.to_excel(writer, sheet_name=scer_fprs_sheet_name, index=False)
        spom_fprs_df.to_excel(writer, sheet_name=spom_fprs_sheet_name, index=False)

        ctrl_sheet_cols = ("Organism", "FPR", "TPR", "Threshold")
        ctrl_sheet.write_row(0, 0, ctrl_sheet_cols)
        ctrl_sheet.write_row(1, 0, ("Calb", "", "", 0))
        ctrl_sheet.write_row(2, 0, ("Scer", "", "", 0))
        ctrl_sheet.write_row(3, 0, ("Spom", "", "", 0))
        ctrl_sheet.add_table(
            0,
            0,
            3,
            3,
            {
                "name": "ControlTable",
                "header_row": True,
                "autofilter": False,
                "columns": [{"header": col} for col in ctrl_sheet_cols],
            },
        )

        calb_fprs_sheet = workbook.get_worksheet_by_name(calb_fprs_sheet_name)
        scer_fprs_sheet = workbook.get_worksheet_by_name(scer_fprs_sheet_name)
        spom_fprs_sheet = workbook.get_worksheet_by_name(spom_fprs_sheet_name)

        calb_fprs_sheet.add_table(
            0,
            0,
            calb_fprs_df.shape[0],
            2,
            {
                "name": "CalbFprTable",
                "header_row": True,
                "autofilter": False,
                "columns": [{"header": col} for col in list(calb_fprs_df.columns.values)],
            },
        )
        scer_fprs_sheet.add_table(
            0,
            0,
            scer_fprs_df.shape[0],
            2,
            {
                "name": "ScerFprTable",
                "header_row": True,
                "autofilter": False,
                "columns": [{"header": col} for col in list(scer_fprs_df.columns.values)],
            },
        )
        spom_fprs_sheet.add_table(
            0,
            0,
            spom_fprs_df.shape[0],
            2,
            {
                "name": "SpomFprTable",
                "header_row": True,
                "autofilter": False,
                "columns": [{"header": col} for col in list(spom_fprs_df.columns.values)],
            },
        )

        ctrl_sheet.write_formula(
            1,
            1,
            "=INDEX(CalbFprTable[FPR], MATCH(ControlTable[[#This Row],[Threshold]], CalbFprTable[Threshold], 0))",
        )
        ctrl_sheet.write_formula(
            1,
            2,
            "=INDEX(CalbFprTable[TPR], MATCH(ControlTable[[#This Row],[Threshold]], CalbFprTable[Threshold], 0))",
        )
        ctrl_sheet.write_formula(
            2,
            1,
            "=INDEX(ScerFprTable[FPR], MATCH(ControlTable[[#This Row],[Threshold]], ScerFprTable[Threshold], 0))",
        )
        ctrl_sheet.write_formula(
            2,
            2,
            "=INDEX(ScerFprTable[TPR], MATCH(ControlTable[[#This Row],[Threshold]], ScerFprTable[Threshold], 0))",
        )
        ctrl_sheet.write_formula(
            3,
            1,
            "=INDEX(SpomFprTable[FPR], MATCH(ControlTable[[#This Row],[Threshold]], SpomFprTable[Threshold], 0))",
        )
        ctrl_sheet.write_formula(
            3,
            2,
            "=INDEX(SpomFprTable[TPR], MATCH(ControlTable[[#This Row],[Threshold]], SpomFprTable[Threshold], 0))",
        )

        orth_sheet = workbook.get_worksheet_by_name(orth_sheet_name)
        orth_sheet.add_table(
            0,
            0,
            orth_df.shape[0],
            orth_df.shape[1] - 1,
            {
                "name": "OrthologsTable",
                "header_row": True,
                "columns": [{"header": col} for col in list(orth_df.columns.values)],
            },
        )
        for org_ix, prefix in enumerate(("Ca", "Sc", "Sp")):
            verdict_col_ix = list(orth_df.columns.values).index("%s verdict" % prefix)
            for row_ix in _xrange(1, orth_df.shape[0] + 1):
                orth_sheet.write_formula(
                    row_ix,
                    verdict_col_ix,
                    '=IF(OrthologsTable[[#This Row],[%s RF G4]]>=\'%s\'!$D$%d, "Yes", "No")'
                    % (prefix, ctrl_sheet_name, org_ix + 2),
                )


# ============================================================
# NEW: pipeline-agnostic "table mode" for glabrata classification
# ============================================================
class TableFeature(object):
    """
    Minimal Feature-like object so COLS_CONFIG can keep working.

    We store:
      - standard_name
      - common_name
      - type
      - description
      - cerevisiae_orthologs (set)
      - feature_qualifier (for no_dubious_features)
      - feature_name (used in some figure name matching)
    """

    __slots__ = (
        "standard_name",
        "common_name",
        "type",
        "description",
        "cerevisiae_orthologs",
        "feature_qualifier",
        "feature_name",
    )

    def __init__(self, standard_name, common_name="", ftype="ORF", description=""):
        self.standard_name = str(standard_name)
        self.common_name = str(common_name) if common_name is not None else ""
        self.type = str(ftype) if ftype else "ORF"
        self.description = str(description) if description is not None else ""
        self.cerevisiae_orthologs = set()
        self.feature_qualifier = ""
        self.feature_name = self.common_name or self.standard_name

    def __len__(self):
        # In table mode, true gene length is stored as record["length"]
        # but COLS_CONFIG uses len(feature). We patch that by storing
        # a proxy in record["feature_len_proxy"] and setting __len__ via closure is not possible.
        # So: keep __len__ as 0, and in table mode we override COLS_CONFIG entry for Length.
        return 0

def _read_any_table(path):
    # supports CSV/TSV; handles "RDF,1" header line in *.feature_table.RDF_1.csv
    with open(path, "r") as fh:
        first = fh.readline().strip()
        second = fh.readline().strip()

    skiprows = 1 if first.startswith("RDF") else 0

    # extension-based default
    if path.lower().endswith((".tsv", ".txt")):
        sep = "\t"
    else:
        # sniff delimiter for "csv" too (some are actually tab-delimited)
        sep = "\t" if ("\t" in second and "," not in second) else ","

    df = pd.read_csv(path, sep=sep, skiprows=skiprows)
    df.columns = [str(c).strip() for c in df.columns]

    # guardrail: catch the classic bad parse early
    if df.shape[1] <= 2 and set(df.columns) <= {"RDF", "1"}:
        raise ValueError(
            f"{path} parsed as columns {list(df.columns)}; looks like RDF header wasn't handled "
            f"or delimiter was wrong."
        )

    return df

def load_bg2_scer_orthology_map(orth_file):
    """
    Expected columns (minimum): cglab_id, scer_id
    Optional: scer_gene_name, Orthogroup, scer_sgdid, etc.
    Returns dict keyed by cglab_id -> dict with scer_id_list + summary fields.
    """
    orth = pd.read_csv(orth_file, sep="\t")
    needed = {"cglab_id", "scer_id"}
    missing = needed - set(orth.columns)
    if missing:
        raise ValueError("Orthology file missing columns: %s" % (sorted(missing),))

    out = {}
    for cglab_id, sub in orth.groupby("cglab_id", dropna=False):
        sub2 = sub.dropna(subset=["scer_id"])
        scer_ids = [str(x) for x in sub2["scer_id"].tolist() if str(x).strip() != ""]
        scer_ids = list(dict.fromkeys(scer_ids))
        first = scer_ids[0] if scer_ids else ""

        first_row = sub2.iloc[0] if len(sub2) else sub.iloc[0]

        out[str(cglab_id)] = {
            "scer_id": first,
            "scer_id_list": scer_ids,
            "scer_gene_name": str(first_row.get("scer_gene_name", "") or ""),
            "Orthogroup": str(first_row.get("Orthogroup", "") or ""),
            "scer_sgdid": str(first_row.get("scer_sgdid", "") or ""),
        }
    return out


def read_glabrata_feature_table_as_records(feature_table_csv, orth_map=None):
    """
    Reads BG2 *.feature_table.RDF_1.csv into "record" dicts compatible with run_pipeline.
    """
    orth_map = orth_map or {}
    df = _read_any_table(feature_table_csv)

    # enforce/normalize columns
    required = {
        "standard_name",
        "coding_length",
        "hits",
        "reads",
        "neighborhood_index",
        "freedom_index",
        "upstream_hits_100",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError("%s missing required columns: %s" % (feature_table_csv, sorted(missing)))

    records = []
    for row in df.itertuples(index=False):
        rd = row._asdict()

        std = str(rd.get("standard_name", ""))
        feat = TableFeature(
            standard_name=std,
            common_name=str(rd.get("common_name", "") or ""),
            ftype=str(rd.get("type", "ORF") or "ORF"),
            description=str(rd.get("description", "") or ""),
        )

        # attach scer orthology if provided
        orth = orth_map.get(std)
        if orth and orth.get("scer_id_list"):
            feat.cerevisiae_orthologs = set(orth.get("scer_id_list", []))

        rec = {"feature": feat}

        # Map length + keep coding_length
        rec["coding_length"] = int(rd.get("coding_length", 0) or 0)
        rec["length"] = rec["coding_length"]

        # Core metrics used by G4:
        rec["hits"] = int(rd.get("hits", 0) or 0)
        rec["reads"] = int(rd.get("reads", 0) or 0)
        rec["neighborhood_index"] = float(rd.get("neighborhood_index", 0.0) or 0.0)
        rec["freedom_index"] = float(rd.get("freedom_index", 0.0) or 0.0)
        rec["upstream_hits_100"] = int(rd.get("upstream_hits_100", 0) or 0)

        # carry-through any extra columns
        for k, v in rd.items():
            if k not in rec:
                rec[k] = v

        # add orthology annotation fields for convenience in outputs
        if orth:
            rec["scer_id"] = orth.get("scer_id", "")
            rec["scer_gene_name"] = orth.get("scer_gene_name", "")
            rec["Orthogroup"] = orth.get("Orthogroup", "")
            rec["scer_sgdid"] = orth.get("scer_sgdid", "")

        # table mode has no ground truth
        rec["ground_truth"] = ""

        records.append(rec)

    return records


def build_scer_training_records_from_dependencies():
    """
    Builds Scer training records exactly like the original script did:
    - Reads Kornmann WildType wig tracks from dependencies
    - analyze_hits against Organisms.cer.feature_db
    - filters ORFs, non-dubious, non-depleted, ignores Organisms.cer.ignored_features
    Returns:
      cer_records_all, training_ess_records, training_non_ess_records
    """
    cer_db = Organisms.cer.feature_db
    all_cer_track_files = glob.glob(Shared.get_dependency("Kornmann", "*WildType*.wig"))
    if not all_cer_track_files:
        raise RuntimeError("No Kornmann *WildType*.wig tracks found in dependencies.")

    all_cer_tracks = [SummaryTable.get_hits_from_wig(fname) for fname in all_cer_track_files]
    cer_wt_combined = sum(all_cer_tracks, [])

    cer_wt_combined_analyzed = list(SummaryTable.analyze_hits(cer_wt_combined, cer_db, 10000).values())
    cer_records = list(only_orfs(no_dubious_features(no_depleted_hit_features(cer_wt_combined_analyzed))))

    ignored = Organisms.cer.ignored_features
    cer_records = [r for r in cer_records if r["feature"].standard_name not in ignored]

    for r in cer_records:
        name = r["feature"].standard_name
        r["ground_truth"] = (
            "Yes"
            if name in Organisms.cer.literature_essentials
            else "No"
            if name in Organisms.cer.literature_non_essentials
            else ""
        )
        r["has_paralog"] = "Yes" if name in Organisms.cer.genes_with_paralogs else ""
        r["filtered_in_training"] = ""

    # Manual curated Scer FP/FN 
    scer_training_fps = set(["S000005081", "S000000510", "S000003278", "S000005250", "S000002448", "S000006163"])
    scer_training_fns = set(
        [
            "S000000370",
            "S000005217",
            "S000003026",
            "S000005194",
            "S000000175",
            "S000002583",
            "S000003092",
            "S000001390",
            "S000003850",
            "S000000223",
            "S000001915",
            "S000004391",
            "S000001295",
            "S000003716",
            "S000003148",
            "S000003073",
            "S000002856",
            "S000003063",
            "S000003440",
            "S000005710",
            "S000000393",
            "S000002248",
            "S000005850",
            "S000000363",
            "S000006169",
            "S000004388",
            "S000006335",
            "S000002653",
            "S000005024",
            "S000001537",
            "S000003293",
            "S000000070",
        ]
    )

    for r in cer_records:
        name = r["feature"].standard_name
        r["filtered_in_training"] = "FP" if name in scer_training_fps else "FN" if name in scer_training_fns else ""

    ess_names = (Organisms.cer.literature_essentials - ignored) - scer_training_fps
    non_names = (Organisms.cer.literature_non_essentials - ignored) - scer_training_fns

    training_ess = [dict(r) for r in cer_records if r["feature"].standard_name in ess_names]
    training_non = [dict(r) for r in cer_records if r["feature"].standard_name in non_names]

    return cer_records, training_ess, training_non


def _table_mode_cols_config():
    """
    In table mode, len(feature) is not meaningful, because TableFeature.__len__ is 0.
    So we clone COLS_CONFIG and override the "Length" column to use record['length'].
    """
    cfg = []
    for c in COLS_CONFIG:
        if c.get("csv_name") == "Length":
            cfg.append({"field_name": "length", "csv_name": "Length", "format": "%d", "local": True})
        else:
            cfg.append(c)
    return cfg


def _table_mode_scer_training_cols_config():
    """
    Scer training sheet (table mode):
      - Standard name = systematic ORF (e.g., YBR089C)  [feature.feature_name]
      - Common name   = gene symbol (e.g., PKC1)       [feature.common_name]
      - SGDID         = SGD identifier (e.g., S000000201) [feature.standard_name]
      - Drop Sc ortholog / Sc std name (not meaningful for Scer training)
    """
    cfg = []
    for c in COLS_CONFIG:
        name = c.get("csv_name", "")

        # drop ortholog columns for Scer training
        if name in ("Sc ortholog", "Sc std name"):
            continue

        # Replace "Standard name" with ORF name (YBR089C)
        if name == "Standard name":
            cfg.append({
                "field_name": "feature",
                "csv_name": "Standard name",
                "format": lambda f: getattr(f, "feature_name", "") or "",
            })
            # add SGDID right after
            cfg.append({
                "field_name": "feature",
                "csv_name": "SGDID",
                "format": lambda f: getattr(f, "standard_name", "") or "",
            })
            continue

        # Replace "Common name" with gene symbol
        if name == "Common name":
            cfg.append({
                "field_name": "feature",
                "csv_name": "Common name",
                "format": lambda f: getattr(f, "common_name", "") or "",
            })
            continue

        # keep everything else
        cfg.append(c)

    return cfg
def main_table_mode_scer_train_gla_classify(args):
    Shared.make_dir(args.out_dir)

    # Build Scer training from dependencies
    _cer_all, train_ess, train_non = build_scer_training_records_from_dependencies()

    # Load orthology map if provided
    orth_map = {}
    if args.orthology_file:
        orth_map = load_bg2_scer_orthology_map(args.orthology_file)

    # Collect glabrata feature tables
    gla_tables = list(args.gla_feature_table or [])
    if args.gla_feature_glob:
        gla_tables.extend(glob.glob(args.gla_feature_glob))
    gla_tables = sorted(set(gla_tables))
    if not gla_tables:
        raise RuntimeError("No glabrata feature tables found. Use --gla-feature-table or --gla-feature-glob")

    # Build datasets
    class_datasets = OrderedDict()
    per_lib_records = OrderedDict()

    for p in gla_tables:
        label = os.path.basename(p)
        label = re.sub(r"\.feature_table\.RDF_1\.csv$", "", label)
        recs = read_glabrata_feature_table_as_records(p, orth_map=orth_map)
        per_lib_records[label] = recs
        class_datasets[label] = (recs,)  # single dataset per label

    # Feature groups (default: original G4)
    feature_groups = OrderedDict(
        (
            ("G4", ("neighborhood_index", "length", "hits", "reads", "freedom_index", "upstream_hits_100")),
        )
    )

    # Classifier factories
    cls_factory = {
        "RF": lambda: RandomForestClassifier(n_estimators=100, random_state=0, n_jobs=int(args.threads)),
        # keep LR option, but off by default unless asked
    }
    if args.include_lr:
        cls_factory["LR"] = lambda: LogisticRegression(max_iter=1000, n_jobs=int(args.threads))

    # Use table-mode cols config
    train_cols_config = _table_mode_scer_training_cols_config()
    class_cols_config = _table_mode_cols_config()

    data = OrderedDict(
        (
            ("Scer training - filtered (table mode)", (train_ess, train_non, class_datasets, OrderedDict())),
        )
    )

    run_pipeline(cls_factory, feature_groups, data,
             train_cols_config,
             class_cols_config,
             args.out_dir, float(args.target_fpr))

    # Optional: combined table (mean score across libraries)
    if args.combine:
        score_key = "RF-G4"
        verdict_key = "RF-G4-verdict"

        rows = []
        for lib, recs in per_lib_records.items():
            for r in recs:
                rows.append(
                    {
                        "standard_name": r["feature"].standard_name,
                        "library": lib,
                        "score": r.get(score_key, np.nan),
                        "verdict": r.get(verdict_key, ""),
                        "hits": r.get("hits", np.nan),
                        "reads": r.get("reads", np.nan),
                        "neighborhood_index": r.get("neighborhood_index", np.nan),
                        "freedom_index": r.get("freedom_index", np.nan),
                        "upstream_hits_100": r.get("upstream_hits_100", np.nan),
                        "length": r.get("length", np.nan),
                        "scer_id": r.get("scer_id", ""),
                        "scer_gene_name": r.get("scer_gene_name", ""),
                        "Orthogroup": r.get("Orthogroup", ""),
                    }
                )
        df = pd.DataFrame(rows)
        comb = (
            df.groupby(["standard_name", "scer_id", "scer_gene_name", "Orthogroup"], dropna=False)
            .agg(score_mean=("score", "mean"), n_libs=("library", "nunique"))
            .reset_index()
        )
        comb.to_csv(os.path.join(args.out_dir, "combined_glabrata_RF-G4.tsv"), sep="\t", index=False)


# ============================================================
# ORIGINAL main() (legacy paper mode) — fully Python3-fixed
# ============================================================
def main():
    # List of folders used for source data:
    # TODO: this should be generalized eventually and refactored with argparse,
    # but currently too much stuff is hardcoded
    output_folder = os.path.join(Shared.get_script_dir(), "output", "predictions")
    alb_hit_file_folder = Shared.get_dependency("albicans", "experiment data", "post evo", "q20m2")
    all_cer_track_files = glob.glob(Shared.get_dependency("Kornmann", "*WildType*.wig"))
    pombe_hit_file = Shared.get_dependency("pombe", "hermes_hits.csv")

    # TODO: should not be hard-coded. Consider using DomainFigures to create the figures
    # as needed, instead of copying them from somewhere else.
    calb_source_figure_folder = "/Users/bermanlab/OneDrive2/OneDrive/Tn Paper/All gene figures/Calb"
    calb_source_figure_list = os.listdir(calb_source_figure_folder)

    scer_source_figure_folder = "/Users/bermanlab/OneDrive2/OneDrive/Tn Paper/All gene figures/Scer"
    scer_source_figure_list = os.listdir(scer_source_figure_folder)

    spom_source_figure_folder = "/Users/bermanlab/OneDrive2/OneDrive/Tn Paper/All gene figures/Spom"
    spom_source_figure_list = os.listdir(spom_source_figure_folder)

    # Set up needed infrastructure
    # The feature DBs:
    alb_db = Organisms.alb.feature_db
    cer_db = Organisms.cer.feature_db
    pom_db = Organisms.pom.feature_db

    # Read cerevisiae hit data:
    all_cer_tracks = [SummaryTable.get_hits_from_wig(fname) for fname in all_cer_track_files]
    cer_wt_combined = sum(all_cer_tracks, [])

    cer_wt_combined_analyzed = list(SummaryTable.analyze_hits(cer_wt_combined, cer_db, 10000).values())
    cer_wt_combined_benchmark_records = list(
        only_orfs(no_dubious_features(no_depleted_hit_features(cer_wt_combined_analyzed)))
    )
    ignored_cer_genes = Organisms.cer.ignored_features
    cer_wt_combined_benchmark_records = [
        r for r in cer_wt_combined_benchmark_records if r["feature"].standard_name not in ignored_cer_genes
    ]

    for record in cer_wt_combined_analyzed:
        rec_name = record["feature"].standard_name
        record["ground_truth"] = (
            "Yes"
            if rec_name in Organisms.cer.literature_essentials
            else "No"
            if rec_name in Organisms.cer.literature_non_essentials
            else ""
        )
        record["has_paralog"] = "Yes" if rec_name in Organisms.cer.genes_with_paralogs else ""

    # Read albicans hit data:
    rdf = 1  # Read depth filter used
    post_hits = SummaryTable.read_hit_files(
        glob.glob(os.path.join(alb_hit_file_folder, "*03*"))
        + glob.glob(os.path.join(alb_hit_file_folder, "*07*"))
        + glob.glob(os.path.join(alb_hit_file_folder, "*11*")),
        rdf,
    )
    post_analyzed = SummaryTable.analyze_hits(sum(post_hits, []), alb_db)

    ignored_alb_genes = Organisms.alb.ignored_features
    ignored_alb_genes |= set(f.standard_name for f in alb_db.get_all_features() if not f.is_orf or "dubious" in f.type.lower())
    for key in list(post_analyzed.keys()):
        if key in ignored_alb_genes:
            del post_analyzed[key]

    SummaryTable.enrich_alb_records(list(post_analyzed.values()))

    # These FPs and FNs were manually curated by us:
    alb_lit_fps = set(
        [
            "C5_02260C_A",
            "C4_07200C_A",
            "C7_01840W_A",
            "C3_07480W_A",
            "C1_09370W_A",
            "C1_06900C_A",
            "C4_04930C_A",
            "C1_14190C_A",
            "C6_00320C_A",
            "C5_05310W_A",
            "C1_11280W_A",
            "C2_05140W_A",
            "C4_05180C_A",
            "C1_14470W_A",
        ]
    )
    alb_lit_fns = set(
        [
            "C4_00610W_A",
            "C2_10210C_A",
            "C1_03600W_A",
            "C7_00700W_A",
            "C4_04730W_A",
            "CR_01740W_A",
            "CR_00710C_A",
            "C2_04760W_A",
            "C1_04380W_A",
            "C4_03440C_A",
            "C1_04090C_A",
            "C4_04180C_A",
            "C1_10860C_A",
            "C5_02900W_A",
            "CR_06420W_A",
            "C2_09370C_A",
            "CR_03240C_A",
            "C1_04330W_A",
            "C7_00890C_A",
            "C2_07100W_A",
            "C2_04220C_A",
            "C6_02840C_A",
            "C5_01720C_A",
            "C1_02230W_A",
            "C1_09870W_A",
            "C2_06540C_A",
            "C1_01790W_A",
            "C3_07550C_A",
            "C5_04600C_A",
            "C3_02960C_A",
            "C5_05190W_A",
            "CR_03430W_A",
            "C1_12510W_A",
            "CR_10140W_A",
            "CR_05620C_A",
            "C1_00700W_A",
            "C4_01190W_A",
            "C1_01490W_A",
            "C4_00130W_A",
            "C7_04230W_A",
            "C1_10210C_A",
            "CR_05030W_A",
            "C1_11400C_A",
            "C1_00060W_A",
            "CR_07580C_A",
            "C7_03940C_A",
            "CR_06640C_A",
            "C4_04090C_A",
            "C4_04850C_A",
            "C1_06230C_A",
            "C1_07970C_A",
            "C2_03180C_A",
        ]
    )

    alb_orths_essential_in_sp_and_sc = SummaryTable.get_calb_ess_in_sc()[0] & SummaryTable.get_calb_ess_in_sp()[0]
    alb_deleted_genes = SummaryTable.get_homann_deletions() | SummaryTable.get_noble_deletions() | SummaryTable.get_sanglard_deletions()

    # These are "systematic", but not manually filtered.
    alb_lit_ess = alb_orths_essential_in_sp_and_sc - alb_deleted_genes
    alb_lit_non_ess = alb_deleted_genes - alb_orths_essential_in_sp_and_sc

    # Used for training:
    filtered_alb_lit_ess = alb_lit_ess - alb_lit_fps
    filtered_alb_lit_non_ess = alb_lit_non_ess - alb_lit_fns

    # TODO: figure out the Roemer thing.
    roemer_ess, roemer_non_ess, omeara_ess, omeara_non_ess = SummaryTable.get_grace_essentials()

    # TODO: generalize for cerevisiae and pombe?
    for record in post_analyzed.values():
        rec_name = record["feature"].standard_name
        record["filtered_in_training"] = "FP" if rec_name in alb_lit_fps else "FN" if rec_name in alb_lit_fns else ""
        record["ground_truth"] = "Yes" if rec_name in alb_lit_ess else "No" if rec_name in alb_lit_non_ess else ""

    # Read the pombe hit data:
    pom_hits = SummaryTable.read_pombe_hit_file(pombe_hit_file, rdf)
    ignored_pom_genes = Organisms.pom.ignored_features
    pom_records = SummaryTable.analyze_hits(pom_hits, pom_db)
    pom_records = {n: f for n, f in pom_records.items() if n not in ignored_pom_genes}
    pom_ess_names = Organisms.pom.literature_essentials
    pom_non_ess_names = Organisms.pom.literature_non_essentials
    pom_ess_names -= ignored_pom_genes
    pom_non_ess_names -= ignored_pom_genes
    pom_ess_records = [pom_records[n] for n in pom_ess_names if n in pom_records]
    pom_non_ess_records = [pom_records[n] for n in pom_non_ess_names if n in pom_records]
    for r in pom_records.values():
        std_name = r["feature"].standard_name
        if std_name in pom_ess_names:
            r["ground_truth"] = "Yes"
        elif std_name in pom_non_ess_names:
            r["ground_truth"] = "No"
        else:
            r["ground_truth"] = ""

        r["has_paralog"] = "Yes" if std_name in Organisms.pom.genes_with_paralogs else ""

    # Set up the classifer:
    cols_config = COLS_CONFIG

    # Define the classification parameters:
    feature_groups = OrderedDict(
        (
            # ("G3", ("neighborhood_index", "length", "hits", "freedom_index", "upstream_hits_100")),
            ("G4", ("neighborhood_index", "length", "hits", "reads", "freedom_index", "upstream_hits_100")),
            # ("G5", ("neighborhood_index", "length", "hits", "reads_ni", "freedom_index", "upstream_hits_100")),
        )
    )

    classifier_factory = {
        # "LR": lambda: LogisticRegression(),
        "RF": lambda: RandomForestClassifier(n_estimators=100, random_state=0),
    }

    pom_db = GenomicFeatures.default_pom_db()
    spom_core_ess = (
        set(
            f.standard_name
            for f in (
                pom_db.get_feature_by_name(r["pombe_ortholog"])
                for r in post_analyzed.values()
                if r["essential_in_pombe"] == r["essential_in_cerevisiae"] == "Yes"
            )
            if f
        )
        - Organisms.pom.ignored_features
    )
    spom_core_training_ess = list(map(dict, (r for n, r in pom_records.items() if n in spom_core_ess)))
    spom_core_non_ess = (
        set(
            f.standard_name
            for f in (
                pom_db.get_feature_by_name(r["pombe_ortholog"])
                for r in post_analyzed.values()
                if r["essential_in_pombe"] == r["essential_in_cerevisiae"] == "No"
            )
            if f
        )
        - Organisms.pom.ignored_features
    )
    spom_core_training_non_ess = list(map(dict, (r for n, r in pom_records.items() if n in spom_core_non_ess)))
    scer_core_ess = (
        set(
            f.standard_name
            for f in (
                cer_db.get_feature_by_name(list(r["feature"].cerevisiae_orthologs)[0])
                for r in post_analyzed.values()
                if r["essential_in_pombe"] == r["essential_in_cerevisiae"] == "Yes"
            )
            if f
        )
        - Organisms.cer.ignored_features
    )
    scer_core_training_ess = list(map(dict, (r for r in cer_wt_combined_analyzed if r["feature"].standard_name in scer_core_ess)))
    scer_core_non_ess = (
        set(
            f.standard_name
            for f in (
                cer_db.get_feature_by_name(list(r["feature"].cerevisiae_orthologs)[0])
                for r in post_analyzed.values()
                if r["essential_in_pombe"] == r["essential_in_cerevisiae"] == "No"
            )
            if f
        )
        - Organisms.cer.ignored_features
    )
    scer_core_training_non_ess = list(map(dict, (r for r in cer_wt_combined_analyzed if r["feature"].standard_name in scer_core_non_ess)))
    calb_core_training_ess = list(map(dict, (r for r in _itervalues(post_analyzed) if r["essential_in_cerevisiae"] == "Yes" and r["essential_in_pombe"] == "Yes")))
    calb_core_training_non_ess = list(map(dict, (r for r in _itervalues(post_analyzed) if r["essential_in_cerevisiae"] == "No" and r["essential_in_pombe"] == "No")))

    copy_calb_class_dataset = lambda: list(map(dict, _itervalues(post_analyzed)))
    copy_scer_com_class_dataset = lambda: list(map(dict, cer_wt_combined_benchmark_records))
    copy_spom_class_dataset = lambda: list(map(dict, (r for k, r in _iteritems(pom_records) if k not in ignored_pom_genes)))

    calb_ortholog_class_dataset = copy_calb_class_dataset()
    scer_ortholog_class_dataset = copy_scer_com_class_dataset()
    spom_ortholog_class_dataset = copy_spom_class_dataset()

    copy_omeara_bench_ess = lambda: [dict(r) for r in post_analyzed.values() if r["feature"].standard_name in omeara_ess]
    copy_omeara_bench_non_ess = lambda: [dict(r) for r in post_analyzed.values() if r["feature"].standard_name in omeara_non_ess]

    copy_calb_bench_ess = lambda: list(map(dict, (r for n, r in _iteritems(post_analyzed) if n in filtered_alb_lit_ess)))
    copy_calb_bench_non_ess = lambda: list(map(dict, (r for n, r in _iteritems(post_analyzed) if n in filtered_alb_lit_non_ess)))

    training_calb_manually_curated_ess = [dict(post_analyzed[n]) for n in filtered_alb_lit_ess if n in post_analyzed]
    training_calb_manually_curated_non_ess = [dict(post_analyzed[n]) for n in filtered_alb_lit_non_ess if n in post_analyzed]

    # These FPs and FNs were manually curated by us:
    scer_training_fps = set(["S000005081", "S000000510", "S000003278", "S000005250", "S000002448", "S000006163"])
    scer_training_fns = set(
        [
            "S000000370",
            "S000005217",
            "S000003026",
            "S000005194",
            "S000000175",
            "S000002583",
            "S000003092",
            "S000001390",
            "S000003850",
            "S000000223",
            "S000001915",
            "S000004391",
            "S000001295",
            "S000003716",
            "S000003148",
            "S000003073",
            "S000002856",
            "S000003063",
            "S000003440",
            "S000005710",
            "S000000393",
            "S000002248",
            "S000005850",
            "S000000363",
            "S000006169",
            "S000004388",
            "S000006335",
            "S000002653",
            "S000005024",
            "S000001537",
            "S000003293",
            "S000000070",
        ]
    )

    for r in cer_wt_combined_benchmark_records:
        name = r["feature"].standard_name
        r["filtered_in_training"] = "FP" if name in scer_training_fps else "FN" if name in scer_training_fns else ""

    # These FPs and FNs were manually curated by us:
    spom_training_fps = set(
        [
            "SPCC330.10",
            "SPBC14C8.14c",
            "SPBC725.17c",
            "SPBC25H2.04c",
            "SPAC144.07c",
            "SPAC31A2.05c",
            "SPBC146.05c",
            "SPAC144.18",
            "SPCC18.12c",
            "SPCC1450.10c",
            "SPBC25H2.06c",
            "SPAC22E12.10c",
            "SPAC4D7.12c",
            "SPBC30D10.02",
            "SPBC3B9.12",
            "SPBC16G5.10",
            "SPAC3G9.12",
            "SPCC16C4.08c",
            "SPAC17A5.13",
            "SPAC1F5.02",
            "SPBC685.05",
            "SPBC21.06c",
            "SPAC806.02c",
            "SPAC16E8.15",
            "SPBC19F8.07",
            "SPCC63.10c",
            "SPBC211.01",
            "SPBC2G2.04c",
            "SPAC31G5.08",
            "SPAC1006.02",
            "SPBC14F5.04c",
            "SPBC9B6.10",
            "SPBC3B9.21",
            "SPAC16A10.06c",
            "SPBC21C3.10c",
            "SPAC56E4.02c",
            "SPAC57A7.11",
            "SPAC144.08",
            "SPBC3B8.01c",
            "SPAC3G9.16c",
            "SPBC16D10.10",
            "SPBC36.12c",
            "SPAC12G12.04",
            "SPAC6F12.05c",
            "SPAC222.03c",
            "SPBC1734.03",
        ]
    )
    spom_training_fns = set(
        [
            "SPAC13A11.03",
            "SPAC1F7.07c",
            "SPCC1442.03",
            "SPBC646.13",
            "SPAC19G12.08",
            "SPAC1F12.07",
            "SPCC191.02c",
            "SPBC405.04c",
            "SPAC890.06",
            "SPAC1B2.03c",
            "SPBC119.05c",
            "SPAP27G11.05c",
            "SPBC23G7.08c",
            "SPBC19F8.03c",
            "SPAC977.17",
            "SPAC30D11.10",
            "SPBC530.13",
            "SPAC18G6.04c",
            "SPBC56F2.10c",
            "SPCC736.06",
            "SPCC622.18",
            "SPCC737.02c",
            "SPBC29A3.01",
            "SPAC6B12.15",
            "SPBC2G2.01c",
            "SPBC119.06",
            "SPBC4B4.03",
            "SPAC328.02",
            "SPAC664.02c",
            "SPAC20H4.07",
            "SPAC8C9.03",
            "SPAC227.01c",
            "SPAC25H1.07",
            "SPAC16E8.13",
            "SPCC594.05c",
            "SPBC776.03",
            "SPAC1952.09c",
            "SPAC2F7.07c",
            "SPAC1F3.10c",
            "SPAC22F3.09c",
            "SPCC830.06",
            "SPCC594.06c",
            "SPAC4G8.10",
            "SPAC4F8.01",
            "SPCC16C4.09",
            "SPAC20G8.10c",
            "SPBC146.13c",
            "SPAC1B3.07c",
            "SPBC3F6.05",
            "SPBC215.05",
            "SPAC17G6.04c",
            "SPBC1706.03",
            "SPBP16F5.07",
            "SPAC644.14c",
            "SPCC61.02",
            "SPAC23C11.08",
            "SPBC1778.06c",
            "SPCC338.14",
            "SPCC18B5.03",
            "SPAC15A10.03c",
            "SPBC32F12.01c",
            "SPBC887.10",
        ]
    )

    for n, r in _iteritems(pom_records):
        r["filtered_in_training"] = "FP" if n in spom_training_fps else "FN" if n in spom_training_fns else ""

    filtered_scer_core_training_ess = [dict(r) for r in scer_core_training_ess if r["feature"].standard_name not in scer_training_fps]
    filtered_scer_core_training_non_ess = [dict(r) for r in scer_core_training_non_ess if r["feature"].standard_name not in scer_training_fns]

    filtered_spom_core_training_ess = [dict(r) for r in spom_core_training_ess if r["feature"].standard_name not in spom_training_fps]
    filtered_spom_core_training_non_ess = [dict(r) for r in spom_core_training_non_ess if r["feature"].standard_name not in spom_training_fns]

    copy_scer_com_bench_ess = lambda: [
        dict(r)
        for r in cer_wt_combined_benchmark_records
        if r["feature"].standard_name in (Organisms.cer.literature_essentials - scer_training_fps)
    ]
    copy_scer_com_bench_non_ess = lambda: [
        dict(r)
        for r in cer_wt_combined_benchmark_records
        if r["feature"].standard_name in (Organisms.cer.literature_non_essentials - scer_training_fns)
    ]

    copy_spom_bench_ess = lambda: list(map(dict, (r for r in pom_ess_records if r["feature"].standard_name not in spom_training_fps)))
    copy_spom_bench_non_ess = lambda: list(map(dict, (r for r in pom_non_ess_records if r["feature"].standard_name not in spom_training_fns)))

    data = OrderedDict(
        (
            (
                "Calb training - manually curated",
                (
                    training_calb_manually_curated_ess,
                    training_calb_manually_curated_non_ess,
                    OrderedDict(
                        (
                            ("Calb-all", (copy_calb_class_dataset(),)),
                            ("Calb-ortholog-only", (calb_ortholog_class_dataset,)),
                        )
                    ),
                    OrderedDict(
                        (
                            ("Scer-com", (copy_scer_com_bench_ess(), copy_scer_com_bench_non_ess())),
                            ("Spom", (copy_spom_bench_ess(), copy_spom_bench_non_ess())),
                            ("Calb", (copy_calb_bench_ess(), copy_calb_bench_non_ess())),
                        )
                    ),
                ),
            ),
            (
                "Scer training - filtered",
                (
                    filtered_scer_core_training_ess,
                    filtered_scer_core_training_non_ess,
                    OrderedDict(
                        (
                            ("Calb", (copy_calb_class_dataset(),)),
                            ("Scer-com", (scer_ortholog_class_dataset,)),
                            ("Spom", (copy_spom_class_dataset(),)),
                        )
                    ),
                    OrderedDict(
                        (
                            ("Scer-com", (copy_scer_com_bench_ess(), copy_scer_com_bench_non_ess())),
                            ("Spom", (copy_spom_bench_ess(), copy_spom_bench_non_ess())),
                            ("Calb", (copy_calb_bench_ess(), copy_calb_bench_non_ess())),
                        )
                    ),
                ),
            ),
            (
                "Spom training - filtered",
                (
                    filtered_spom_core_training_ess,
                    filtered_spom_core_training_non_ess,
                    OrderedDict(
                        (
                            ("Calb", (copy_calb_class_dataset(),)),
                            ("Scer-com", (copy_scer_com_class_dataset(),)),
                            ("Spom", (spom_ortholog_class_dataset,)),
                        )
                    ),
                    OrderedDict(
                        (
                            ("Scer-com", (copy_scer_com_bench_ess(), copy_scer_com_bench_non_ess())),
                            ("Spom", (copy_spom_bench_ess(), copy_spom_bench_non_ess())),
                            ("Calb", (copy_calb_bench_ess(), copy_calb_bench_non_ess())),
                        )
                    ),
                ),
            ),
        )
    )

    run_pipeline(classifier_factory, feature_groups, data, cols_config, cols_config, output_folder, 0.1)

    # ---- everything below kept, Python3-fixed ----
    inter_col_config = list(cols_config)
    inter_col_config[-2:-2] = [
        {"field_name": "ground_truth", "csv_name": "Ess. ground truth"},
        {"field_name": "RF-G4", "csv_name": "RF - G4", "format": "%.3f"},
        {"field_name": "RF-G4-verdict", "csv_name": "RF - G4 - ess. verdict"},
    ]

    score_to_use = "RF-G4"
    verdict_to_use = "%s-verdict" % score_to_use
    train_score_to_use = "train-%s" % score_to_use

    # Draw Venns:
    calb_records = list({r["feature"].standard_name: r for r in calb_ortholog_class_dataset}.items())
    calb_records_only = [p[1] for p in calb_records]
    orthologs = set(k for k, r in post_analyzed.items() if r["pombe_ortholog"] and r["feature"].cerevisiae_orthologs)
    draw_venn(
        "Essential orthologs",
        os.path.join(output_folder, "orthologs_calb_ess.png"),
        [
            set(k for k, r in calb_records if r[verdict_to_use] == "Yes") & orthologs,
            set(k for k, r in calb_records if r["essential_in_cerevisiae"] == "Yes") & orthologs,
            set(k for k, r in calb_records if r["essential_in_pombe"] == "Yes") & orthologs,
        ],
        ["Calb", "Scer", "Spom"],
    )

    draw_venn(
        "Non-Essential orthologs",
        os.path.join(output_folder, "orthologs_calb_non_ess.png"),
        [
            set(k for k, r in calb_records if r[verdict_to_use] == "No") & orthologs,
            set(k for k, r in calb_records if r["essential_in_cerevisiae"] == "No") & orthologs,
            set(k for k, r in calb_records if r["essential_in_pombe"] == "No") & orthologs,
        ],
        ["Calb", "Scer", "Spom"],
    )

    calb_essentials = set(k for k, r in post_analyzed.items() if r["essential_in_albicans_grace_omeara"] and r["essential_in_mitchell"])
    draw_venn(
        "Essential in Calb",
        os.path.join(output_folder, "tn_omeara_mitchell_calb_ess.png"),
        [
            set(k for k, r in calb_records if r[verdict_to_use] == "Yes") & calb_essentials,
            set(k for k, r in calb_records if r["essential_in_albicans_grace_omeara"] == "Yes") & calb_essentials,
            set(k for k, r in calb_records if r["essential_in_mitchell"] == "Yes") & calb_essentials,
        ],
        ["Tn", "O'Meara", "Mitchell"],
    )

    draw_venn(
        "Cmp. with Mitchell ess.",
        os.path.join(output_folder, "tn_mitchell_calb_ess.png"),
        [
            set(k for k, r in calb_records if r[verdict_to_use] == "Yes" and r["essential_in_mitchell"]),
            set(k for k, r in calb_records if r["essential_in_mitchell"] == "Yes"),
        ],
        ["Tn", "Mitchell"],
    )

    draw_venn(
        "Cmp. with O'Meara ess.",
        os.path.join(output_folder, "tn_omeara_calb_ess.png"),
        [
            set(k for k, r in calb_records if r[verdict_to_use] == "Yes" and r["essential_in_albicans_grace_omeara"]),
            set(k for k, r in calb_records if r["essential_in_albicans_grace_omeara"] == "Yes"),
        ],
        ["Tn", "O'Meara"],
    )

    draw_venn(
        "Cmp. with O'Meara non-ess.",
        os.path.join(output_folder, "tn_omeara_calb_non_ess.png"),
        [
            set(k for k, r in calb_records if r[verdict_to_use] == "No" and r["essential_in_albicans_grace_omeara"]),
            set(k for k, r in calb_records if r["essential_in_albicans_grace_omeara"] == "No"),
        ],
        ["Tn", "O'Meara"],
    )

    def dump_records(records, source_figure_folder, source_figure_list, inter_out_folder, add_score=False):
        target_figure_folder = os.path.join(inter_out_folder, "figures")
        if not os.path.exists(target_figure_folder):
            Shared.make_dir(target_figure_folder)

        for record in records:
            feature = record["feature"]
            gene_name = getattr(feature, "feature_name", feature.standard_name)
            for source_figure in source_figure_list:
                if gene_name in source_figure:
                    break
            else:
                gene_name = None

            # TODO: the score is hard-coded :(
            dest = (
                target_figure_folder
                if not add_score
                else os.path.join(target_figure_folder, "%.3f_%s" % (record["train-" + score_to_use], source_figure))
            )
            shutil.copy(os.path.join(source_figure_folder, source_figure), dest)

        with pd.ExcelWriter(os.path.join(inter_out_folder, "genes.xlsx")) as excel_writer:
            sheet = SummaryTable.write_data_to_data_frame(records, inter_col_config)
            sheet.to_excel(excel_writer, sheet_name="Main", index=False)

    for ca_status, sc_status, sp_status in product(("Yes", "No"), repeat=3):
        inter_out_folder = os.path.join(output_folder, "Ca-%s Sc-%s Sp-%s" % (ca_status, sc_status, sp_status))
        if not os.path.exists(inter_out_folder):
            os.mkdir(inter_out_folder)

        records = [
            r
            for r in calb_records_only
            if r["essential_in_cerevisiae"] == sc_status and r["essential_in_pombe"] == sp_status and r[verdict_to_use] == ca_status
        ]

        dump_records(records, calb_source_figure_folder, calb_source_figure_list, inter_out_folder)

    for tn_status, omeara_status in product(("Yes", "No"), repeat=2):
        inter_out_folder = os.path.join(output_folder, "Mitchell-Yes Calb-Tn-%s OMeara-%s" % (tn_status, omeara_status))
        if not os.path.exists(inter_out_folder):
            os.mkdir(inter_out_folder)

        records = [
            r
            for r in calb_records_only
            if r["essential_in_mitchell"] == "Yes"
            and r["essential_in_albicans_grace_omeara"] == omeara_status
            and r[verdict_to_use] == tn_status
        ]

        dump_records(records, calb_source_figure_folder, calb_source_figure_list, inter_out_folder)

    inter_col_config[-2:-2] = [
        {"field_name": "train-RF-G4", "csv_name": "RF - G4", "format": "%.3f"},
        {"field_name": "train-RF-G4-verdict", "csv_name": "RF - G4 - ess. verdict"},
    ]

    # Print out the conflicting annotations in the benchmark datasets:
    for training_label, dataset in _iteritems(data):
        training_ess, training_non_ess = dataset[:2]

        dump_records(
            training_ess,
            {"Calb": calb_source_figure_folder, "Scer": scer_source_figure_folder, "Spom": spom_source_figure_folder}[training_label[:4]],
            {"Calb": calb_source_figure_list, "Scer": scer_source_figure_list, "Spom": spom_source_figure_list}[training_label[:4]],
            os.path.join(output_folder, "classification - %s" % training_label, "Training essentials"),
            add_score=True,
        )
        dump_records(
            training_non_ess,
            {"Calb": calb_source_figure_folder, "Scer": scer_source_figure_folder, "Spom": spom_source_figure_folder}[training_label[:4]],
            {"Calb": calb_source_figure_list, "Scer": scer_source_figure_list, "Spom": spom_source_figure_list}[training_label[:4]],
            os.path.join(output_folder, "classification - %s" % training_label, "Training non-essentials"),
            add_score=True,
        )

    # Ortholog megatable:
    index_dataset = lambda dataset: {r["feature"].standard_name: r for r in dataset}
    ortholog_records = []
    indexed_calb = index_dataset(calb_ortholog_class_dataset)
    indexed_scer = index_dataset(scer_ortholog_class_dataset)
    indexed_spom = index_dataset(spom_ortholog_class_dataset)

    # The left-join implementation
    orthologs = [Organisms.get_orths_by_name(n) for n in indexed_calb.keys()]
    sc_first = list(indexed_scer.values())[0]
    sp_first = list(indexed_spom.values())[0]
    sc_rec_items_to_add = list(zip(list(sc_first.keys()), repeat("")))
    sp_rec_items_to_add = list(zip(list(sp_first.keys()), repeat("")))

    for (calb_orth, scer_orth, spom_orth) in orthologs:
        calb_record = indexed_calb.get(calb_orth.standard_name)
        scer_record = indexed_scer.get(scer_orth.standard_name if scer_orth else "")
        spom_record = indexed_spom.get(spom_orth.standard_name if spom_orth else "")

        combined_record = {}
        for prefix, record in zip(("ca", "sc", "sp"), (calb_record, scer_record, spom_record)):
            if record is not None:
                items_to_add = list(record.items())
            else:
                if prefix == "sc":
                    items_to_add = sc_rec_items_to_add
                elif prefix == "sp":
                    items_to_add = sp_rec_items_to_add
                else:
                    items_to_add = []
            for key, value in items_to_add:
                combined_record["%s_%s" % (prefix, key)] = value

        ortholog_records.append(combined_record)

    orth_cols_config = ORTH_COLS_CONFIG

    def get_roc_curve_from_data(data_dict, name):
        training_ess = data_dict[name][0]
        training_non_ess = data_dict[name][1]
        return roc_curve(
            [1] * len(training_ess) + [0] * len(training_non_ess),
            [r[train_score_to_use] for r in training_ess] + [r[train_score_to_use] for r in training_non_ess],
        )

    write_ortholog_excel(
        SummaryTable.write_data_to_data_frame(ortholog_records, orth_cols_config),
        pd.DataFrame(
            dict(zip(("FPR", "TPR", "Threshold"), get_roc_curve_from_data(data, "Calb training - manually curated"))),
            columns=("FPR", "TPR", "Threshold"),
        ),
        pd.DataFrame(
            dict(zip(("FPR", "TPR", "Threshold"), get_roc_curve_from_data(data, "Scer training - filtered"))),
            columns=("FPR", "TPR", "Threshold"),
        ),
        pd.DataFrame(
            dict(zip(("FPR", "TPR", "Threshold"), get_roc_curve_from_data(data, "Spom training - filtered"))),
            columns=("FPR", "TPR", "Threshold"),
        ),
        os.path.join(output_folder, "orthologs_combined_ex.xlsx"),
    )


# ============================================================
# Original configs (unchanged)
# ============================================================
COLS_CONFIG = [
    {"field_name": "feature", "csv_name": "Standard name", "format": lambda f: f.standard_name},
    {"field_name": "feature", "csv_name": "Common name", "format": lambda f: f.common_name},
    {
        "field_name": "feature",
        "csv_name": "Sc ortholog",
        "format": lambda f: ",".join(f.cerevisiae_orthologs) if hasattr(f, "cerevisiae_orthologs") else "",
    },
    {
        "field_name": "feature",
        "csv_name": "Sc std name",
        "format": lambda f: Organisms.cer.feature_db.get_feature_by_name(list(f.cerevisiae_orthologs)[0]).standard_name
        if len(getattr(f, "cerevisiae_orthologs", [])) > 0
        else "",
    },
    {"field_name": "cer_fitness", "csv_name": "Sc fitness"},
    {"field_name": "essential_in_cerevisiae", "csv_name": "Essential in Sc"},
    {"field_name": "cer_synthetic_lethal", "csv_name": "SL in Sc"},
    {"field_name": "pombe_ortholog", "csv_name": "Sp ortholog"},
    {"field_name": "essential_in_pombe", "csv_name": "Essential in Sp"},
    {"field_name": "essential_in_albicans_grace_omeara", "csv_name": "Essential in O'Meara"},
    {"field_name": "essential_in_mitchell", "csv_name": "Essential in Mitchell"},
    {"field_name": "deleted_in_calb", "csv_name": "Deleted in Calb"},
    {"field_name": "filtered_in_training", "csv_name": "Filtered in training"},
    {"field_name": "feature", "csv_name": "Length", "format": lambda f: len(f)},
    {"field_name": "hits", "csv_name": "Hits", "format": "%d", "local": True},
    {"field_name": "reads", "csv_name": "Reads", "format": "%d", "local": True},
    {"field_name": "max_free_region", "csv_name": "Longest free interval", "format": "%d", "local": True},
    {"field_name": "freedom_index", "csv_name": "Freedom index", "format": "%.2f", "local": True},
    {"field_name": "neighborhood_index", "csv_name": "Neighborhood index", "format": "%.3f", "local": True},
    {"field_name": "upstream_hits_100", "csv_name": "Upstream hits 100", "format": "%d", "local": True},
    {"field_name": "upstream_hits_50", "csv_name": "Upstream hits 50", "format": "%d", "local": True},
    {"field_name": "n_term_hits_100", "csv_name": "N-term hits 100", "format": "%d", "local": True},
    {"field_name": "feature", "csv_name": "Type", "format": lambda f: f.type},
    {"field_name": "feature", "csv_name": "Description", "format": lambda f: f.description},
]

ORTH_COLS_CONFIG = [
    {"field_name": "ca_feature", "csv_name": "Ca standard name", "format": lambda f: f.standard_name},
    {"field_name": "ca_feature", "csv_name": "Ca common name", "format": lambda f: f.common_name},
    {"field_name": "sc_feature", "csv_name": "Sc standard name", "format": lambda f: f.standard_name if f else ""},
    {"field_name": "sc_feature", "csv_name": "Sc common name", "format": lambda f: (f.common_name if f else "") or (f.feature_name if f else "")},
    {"field_name": "sp_feature", "csv_name": "Sp standard name", "format": lambda f: f.standard_name if f else ""},
    {"field_name": "sp_feature", "csv_name": "Sp common name", "format": lambda f: f.common_name if f else ""},
    {"field_name": "ca_essential_in_albicans_grace_omeara", "csv_name": "Essential in O'Meara"},
    {"field_name": "ca_essential_in_mitchell", "csv_name": "Essential in Mitchell"},
    {"field_name": "ca_deleted_in_calb", "csv_name": "Deleted in Calb"},
    {"field_name": "ca_has_paralog", "csv_name": "Paralog in Calb"},
    {"field_name": "ca_essential_in_cerevisiae", "csv_name": "Essential in Sc"},
    {"field_name": "ca_cer_synthetic_lethal", "csv_name": "SL in Sc"},
    {"field_name": "sc_has_paralog", "csv_name": "Paralog in Sc"},
    {"field_name": "ca_essential_in_pombe", "csv_name": "Essential in Sp"},
    {"field_name": "sp_has_paralog", "csv_name": "Paralog in Spom"},
    {"field_name": "ca_feature", "csv_name": "Ca length", "format": lambda f: len(f)},
    {"field_name": "sc_feature", "csv_name": "Sc length", "format": lambda f: len(f) if f else ""},
    {"field_name": "sp_feature", "csv_name": "Sp length", "format": lambda f: len(f) if f else ""},
    {"field_name": "ca_hits", "csv_name": "Ca hits", "format": "%d"},
    {"field_name": "sc_hits", "csv_name": "Sc hits", "format": "%d"},
    {"field_name": "sp_hits", "csv_name": "Sp hits", "format": "%d"},
    {"field_name": "ca_reads", "csv_name": "Ca reads", "format": "%d"},
    {"field_name": "sc_reads", "csv_name": "Sc reads", "format": "%d"},
    {"field_name": "sp_reads", "csv_name": "Sp reads", "format": "%d"},
    {"field_name": "ca_max_free_region", "csv_name": "Ca longest free interval", "format": "%d"},
    {"field_name": "sc_max_free_region", "csv_name": "Sc longest free interval", "format": "%d"},
    {"field_name": "sp_max_free_region", "csv_name": "Sp longest free interval", "format": "%d"},
    {"field_name": "ca_freedom_index", "csv_name": "Ca freedom index", "format": "%.2f"},
    {"field_name": "sc_freedom_index", "csv_name": "Sc freedom index", "format": "%.2f"},
    {"field_name": "sp_freedom_index", "csv_name": "Sp freedom index", "format": "%.2f"},
    {"field_name": "ca_neighborhood_index", "csv_name": "Ca neighborhood index", "format": "%.3f"},
    {"field_name": "sc_neighborhood_index", "csv_name": "Sc neighborhood index", "format": "%.3f"},
    {"field_name": "sp_neighborhood_index", "csv_name": "Pp neighborhood index", "format": "%.3f"},
    {"field_name": "ca_upstream_hits_100", "csv_name": "Ca upstream hits 100", "format": "%d"},
    {"field_name": "sc_upstream_hits_100", "csv_name": "Sc upstream hits 100", "format": "%d"},
    {"field_name": "sp_upstream_hits_100", "csv_name": "Sp upstream hits 100", "format": "%d"},
    {"field_name": "ca_RF-G4", "csv_name": "Ca RF G4", "format": "%.3f"},
    {"field_name": "sc_RF-G4", "csv_name": "Sc RF G4", "format": "%.3f"},
    {"field_name": "sp_RF-G4", "csv_name": "Sp RF G4", "format": "%.3f"},
    {"field_name": "ca_feature", "csv_name": "Ca description", "format": lambda f: f.description},
    {"field_name": "sc_feature", "csv_name": "Sc description", "format": lambda f: f.description if f else ""},
    {"field_name": "sp_feature", "csv_name": "Sp description", "format": lambda f: f.description if f else ""},
]


# ============================================================
# CLI wrapper (new; does NOT remove legacy main)
# ============================================================
def _parse_args():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--mode",
        default="scer-train-gla-classify",
        choices=["scer-train-gla-classify", "legacy-paper"],
        help="Run new table mode or the original hard-coded legacy main().",
    )
    ap.add_argument("--out-dir", default=None, help="Output directory (table mode).")
    ap.add_argument("--threads", type=int, default=1, help="Threads (RF n_jobs).")
    ap.add_argument("--target-fpr", type=float, default=0.10, help="Target FPR used for verdict threshold selection.")
    ap.add_argument("--include-lr", action="store_true", help="Also run LogisticRegression in table mode.")

    ap.add_argument("--gla-feature-table", action="append", default=[], help="Path to one glabrata feature_table.RDF_1.csv (repeatable).")
    ap.add_argument("--gla-feature-glob", default=None, help="Glob for glabrata feature_table.RDF_1.csv files.")
    ap.add_argument("--orthology-file", default=None, help="BG2↔Scer orthology TSV (e.g., C_glabrata_BG2_S_cerevisiae_orthologs.txt).")
    ap.add_argument("--combine", action="store_true", help="Also write combined_glabrata_RF-G4.tsv (mean score across libraries).")
    ap.add_argument("--overwrite", action="store_true", help="Allow overwriting an existing output directory / outputs.",)

    return ap.parse_args()

def prepare_out_dir(out_dir: str, overwrite: bool):
    import os
    import shutil

    if os.path.exists(out_dir):
        if not overwrite:
            raise SystemExit(
                f"ERROR: output directory already exists:\n  {out_dir}\n"
                f"Re-run with --overwrite to replace it."
            )
        shutil.rmtree(out_dir)

    os.makedirs(out_dir, exist_ok=True)

if __name__ == "__main__":
    args = _parse_args()

    if args.mode == "legacy-paper":
        main()
    else:
        # ---- table mode sanity checks ----
        if not args.out_dir:
            raise SystemExit("--out-dir is required for table mode.")

        if not args.gla_feature_table and not args.gla_feature_glob:
            raise SystemExit(
                "ERROR: provide at least one --gla-feature-table or --gla-feature-glob"
            )

        prepare_out_dir(args.out_dir, args.overwrite)

        main_table_mode_scer_train_gla_classify(args)
