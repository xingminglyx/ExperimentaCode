import pandas as pd
import numpy as np
from sklearn import metrics
import json

def run_test():
  results = pd.read_csv("features/nci1.csv", header=0, index_col=0)

  confusion = pd.DataFrame(np.zeros((2, 2)), columns=("LUAD", "LUSC"), index=("LUAD", "LUSC"))
  y_true = []
  y_pred = []

  with open("./patients_split_list.json") as f:
    patients_split_list = json.load(f)["testing"]
    for disease in ("LUAD", "LUSC"):
      for dataset in patients_split_list[disease].keys():
        for patient in patients_split_list[disease][dataset]:
          # dis1 = np.sum(np.abs(results.loc["LUAD_template", ] - results.loc[patient, ]))
          # dis2 = np.sum(np.abs(results.loc["LUSC_template", ] - results.loc[patient, ]))
          dis1 = np.linalg.norm(results.loc["LUAD_template", ] - results.loc[patient, ])
          dis2 = np.linalg.norm(results.loc["LUSC_template", ] - results.loc[patient, ])
          y_true.append(disease)
          if dis1 < dis2:
            y_pred.append("LUAD")
            if disease == "LUAD":
              confusion.loc["LUAD", "LUAD"] += 1
            else:
              confusion.loc["LUAD", "LUSC"] += 1
          elif dis1 > dis2:
            y_pred.append("LUSC")
            if disease == "LUSC":
              confusion.loc["LUSC", "LUSC"] += 1
            else:
              confusion.loc["LUSC", "LUAD"] += 1
          else:
            print(f"{patient} unknow!")

  print(confusion)
  evaluation(y_true, y_pred)

def avg_score(scores):
  std = np.std(scores)
  mean = np.mean(scores)
  return f"{round(mean, 4)} ± {round(1.96 * std / np.sqrt(len(scores)), 4)}"

def class_from_prob(x, operating_point=0.5):
  x[x >= operating_point] = 1
  x[x < operating_point] = 0
  return x

def evaluation(y_true, y_pred, y_prob):
  accuracy = metrics.accuracy_score(y_true, y_pred)
  precision = metrics.precision_score(y_true, y_pred, pos_label=1)
  recall = metrics.recall_score(y_true, y_pred, pos_label=1)
  F1_score = metrics.f1_score(y_true, y_pred, pos_label=1)
  MCC = metrics.matthews_corrcoef(y_true, y_pred)
  kappa = metrics.cohen_kappa_score(y_true, y_pred)
  auc = metrics.roc_auc_score(y_true, y_prob)
  # fpr, tpr, thresholds = metrics.roc_curve(y_true, y_prob, pos_label=1)
  # auc = metrics.auc(fpr, tpr)
  return accuracy, precision, recall, F1_score, MCC, kappa, auc

def multi_evaluation(y_true, y_pred, y_prob):
  accuracy = metrics.accuracy_score(y_true, y_pred)
  F1_score_micro = metrics.f1_score(y_true, y_pred, average='micro')
  F1_score_macro = metrics.f1_score(y_true, y_pred, average='macro')
  MCC = metrics.matthews_corrcoef(y_true, y_pred)
  auc = metrics.roc_auc_score(y_true, y_prob, multi_class='ovo')
  return accuracy, F1_score_micro, F1_score_macro, MCC, auc

def add_results(all_results, method, y_true, y_pred, y_prob, multi = False):
  if multi:
    accuracy, F1_score_micro, F1_score_macro, MCC, auc = multi_evaluation(y_true, y_pred, y_prob)
    all_results[method]["accuracy"].append(accuracy)
    all_results[method]["F1_micro"].append(F1_score_micro)
    all_results[method]["F1_macro"].append(F1_score_macro)
    all_results[method]["MCC"].append(MCC)
    all_results[method]["auc"].append(auc)
    return all_results
  accuracy, precision, recall, F1, mcc, kappa, auc = evaluation(y_true, y_pred, y_prob)
  all_results[method]["accuracy"].append(accuracy)
  all_results[method]["precision"].append(precision)
  all_results[method]["recall"].append(recall)
  all_results[method]["F1"].append(F1)
  all_results[method]["mcc"].append(mcc)
  all_results[method]["kappa"].append(kappa)
  all_results[method]["auc"].append(auc)
  return all_results
