/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.mines.processing.experiments;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import org.mines.processing.utils.ResultsTable;
import org.mines.processing.utils.Result;
import org.mines.processing.utils.ComputeResult;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import org.data2semantics.mustard.kernels.Kernel;
import org.data2semantics.mustard.kernels.KernelUtils;
import org.data2semantics.mustard.kernels.SparseVector;
import org.data2semantics.mustard.kernels.data.RDFData;
import org.data2semantics.mustard.learners.Prediction;
import org.data2semantics.mustard.learners.evaluation.AUCROC;
import org.data2semantics.mustard.learners.evaluation.Accuracy;
import org.data2semantics.mustard.learners.evaluation.EvaluationFunction;
import org.data2semantics.mustard.learners.evaluation.F1;
import org.data2semantics.mustard.learners.evaluation.utils.EvaluationUtils;
import org.data2semantics.mustard.learners.libsvm.LibSVM;
import org.data2semantics.mustard.learners.libsvm.LibSVMParameters;
import org.data2semantics.mustard.rdf.DataSetUtils;
import org.data2semantics.mustard.rdf.RDFDataSet;
import org.data2semantics.mustard.rdf.RDFFileDataSet;
import org.openrdf.model.Resource;
import org.openrdf.model.Statement;
import org.openrdf.model.Value;
import org.openrdf.rio.RDFFormat;
import org.data2semantics.mustard.learners.libsvm.LibSVMModel;

/**
 *
 * @author yassin
 */
public class Classification {

    public static void main(String[] args) throws ClassNotFoundException, InstantiationException, IllegalAccessException, IllegalArgumentException, InvocationTargetException, NoSuchMethodException, FileNotFoundException {

        PrintStream out = new PrintStream(new FileOutputStream("output-classification.txt", false));
        System.setOut(out);
        System.setErr(out);
        
        //load rdf graph in memory
        System.out.println("loading graph in memory...");
        long loadstart = System.currentTimeMillis();
        String datasetPath = args[0];
        //String datasetPath = "../preprocessing/outputs/directed_subgraph";
        RDFDataSet tripleStore = new RDFFileDataSet(datasetPath, RDFFormat.RDFXML);
        long loadend = System.currentTimeMillis();
        System.out.println("graph loaded in: "+ (loadend - loadstart)+" ms");

        ResultsTable resTable = new ResultsTable();
        resTable.setDigits(3);
        resTable.setSignificanceTest(ResultsTable.SigTest.PAIRED_TTEST);
        resTable.setpValue(0.05);
        resTable.setShowStdDev(true);
        resTable.setLatex(true);

        // general blacklist
        List<Resource> allinstances = new ArrayList<>();
        List<Value> alllabels = new ArrayList<>();
        for (int num = 0; num <= 9; num++) {
            for (Statement stmt : tripleStore.getStatementsFromStrings(null, "http://bio2rdf.org/pharmgkb_vocabulary:train_subset_" + num, null)) {
                allinstances.add(stmt.getSubject());
                alllabels.add(stmt.getObject());
            }
        }
        for (Statement stmt : tripleStore.getStatementsFromStrings(null, "http://bio2rdf.org/pharmgkb_vocabulary:train_subset_final", null)) {
            allinstances.add(stmt.getSubject());
            alllabels.add(stmt.getObject());
        }
        for (Statement stmt : tripleStore.getStatementsFromStrings(null, "http://bio2rdf.org/pharmgkb_vocabulary:test_set", null)) {
            allinstances.add(stmt.getSubject());
            alllabels.add(stmt.getObject());
        }
        List<Statement> blackList = DataSetUtils.createBlacklist(tripleStore, allinstances, alllabels);
        
        

        // run performance estimation over the 10 random training sets
        for (int num = 10; num <= 10; num++) {
            if (num != 10) {
                resTable.newRow("sample " + num);
            } else {
                resTable.newRow("final model");
            }

            List<Resource> instances = new ArrayList<>();
            List<Value> labels = new ArrayList<>();
            String trainURI = (num != 10) ? "http://bio2rdf.org/pharmgkb_vocabulary:train_subset_" + num : "http://bio2rdf.org/pharmgkb_vocabulary:train_subset_final";
            long tic = System.currentTimeMillis();
            List<Statement> trainStmts = new ArrayList<>();
            trainStmts = tripleStore.getStatementsFromStrings(null, trainURI, null);
            for (Statement stmt : trainStmts) {
                instances.add(stmt.getSubject());
                labels.add(stmt.getObject());
            }
            double[] targets = LibSVM.createTargets(labels);
            RDFData data = new RDFData(tripleStore, instances, blackList);

            /// kernels
            String[] kernels = {"RDFGraphListWLSubTreeKernel", "RDFTreeWLSubTreeKernel", "RDFGraphListWalkCountKernel", "RDFTreeWalkCountKernel",
                "RDFGraphListWalkCountApproxKernelMkII", "RDFTreeWalkCountIDEQApproxKernelMkII", "RDFRootWalkCountKernel", "RDFRootWLSubTreeKernel",
                "RDFGraphListWLSubTreeApproxKernel", "RDFTreeWLSubTreeIDEQApproxKernel", "RDFRootWalkCountIDEQApproxKernel", "RDFRootWLSubTreeIDEQApproxKernel"};
             
            // kernels params
            // For details about the last 3 , see https://github.com/Data2Semantics/mustard/blob/master/mustard-kernels/src/main/java/org/data2semantics/mustard/kernels/graphkernels/graphlist/WLSubTreeApproxKernel.java
            int[] pathLength_values = {0, 4, 6, 8, 10, 12, 15};
            int[] depth_values = {4, 6, 8, 10, 12, 15};
            boolean[] inference_values = {false};
            boolean[] normalize_values = {true};
            int[] minFreqs = {0,2,4,8,16,32};
            int[] maxPrevNBHs = {0, 1, 5};
            int[] maxLabelCards = {1, 10000};
            // svm params
            double[] cs = {1, 10, 100, 1000};
            //double[] nus = {0.1, 0.3, 0.5, 0.7, 0.9};
            LibSVMParameters svmParms = new LibSVMParameters(LibSVMParameters.C_SVC, cs);
            svmParms.setNumFolds(10);
            svmParms.setEvalFunction(new F1());
            svmParms.setVerbosity(LibSVMParameters.VERBOSITY_NONE);
            
            if (num != 10) {
                System.out.println("Computing kernels for sample " + num + "...");
            } else {
                System.out.println("Computing kernels for final model...");
            }

            Map<Kernel, double[][]> currentKernels = new HashMap<>();
            for (String kernelName : kernels) {
                List<Object[]> paramsCombinations = createKernelParams(kernelName, pathLength_values, depth_values, inference_values, normalize_values, minFreqs, maxPrevNBHs, maxLabelCards);
                for (Object[] params : paramsCombinations) {
                    ComputeResult result = compute(kernelName, params, data);
                    currentKernels.put(result.getKernel(), result.getMatrix());
                }
            }
            
            List<Double> targetsAsList = new ArrayList<>();
                for (int i = 0; i < targets.length; i++) {
                    targetsAsList.add(targets[i]);
                }

            // cross-validate
            // init results
            List<EvaluationFunction> evalFuncs = new ArrayList<EvaluationFunction>();
            evalFuncs.add(new Accuracy());
            evalFuncs.add(new F1());
            evalFuncs.add(new AUCROC());
            long[] seeds = {11, 21, 31, 41, 51, 61, 71, 81, 91, 101};
            List<Result> results = new ArrayList<>();
            for (EvaluationFunction evalFunc : evalFuncs) {
                Result res = new Result(evalFunc);
                double[] resA = new double[seeds.length]; // add a new empty array with the length of the amount of seeds (i.e. the number of repetitions of the experiment).
                res.setScores(resA);
                results.add(res);
            }
            Result compR = new Result();
            compR.setLabel("Computation time");
            results.add(compR);
            if (num != 10) {
                System.out.println("Performing cross-validation for sample " + num + "...");
            } else {
                System.out.println("Performing cross-validation for final model...");
            }

            for (int j = 0; j < seeds.length; j++) {
                for (Kernel k : currentKernels.keySet()) {
                    currentKernels.put(k, KernelUtils.shuffle(currentKernels.get(k), seeds[j]));
                }
                
                Collections.shuffle(targetsAsList, new Random(seeds[j]));

                double[] newtarget = new double[targetsAsList.size()];
                for (int i = 0; i < newtarget.length; i++) {
                    newtarget[i] = targetsAsList.get(i);
                }
                Prediction[] pred = LibSVM.crossValidateWithMultipleKernels(currentKernels, newtarget, svmParms, svmParms.getNumFolds());
                for (Result res : results) {
                    if (res.getEval() != null) {
                        res.getScores()[j] = res.getEval().computeScore(newtarget, pred);
                    }
                }

            }
            long toc = System.currentTimeMillis();
            double[] comp = {toc - tic};
            compR.setScores(comp);
            for (Result res : results) {
                resTable.addResult(res);
            }

        }

        resTable.addCompResults(resTable.getBestResults());
        System.out.println(resTable);

        for (int i = 0; i <= 0; i++) {
            List<Resource> instances = new ArrayList<>();
            List<Value> labels = new ArrayList<>();
            System.out.println("training final model for prediction (using training set and test set)...");
            long tic = System.currentTimeMillis();

            String trainURI = "http://bio2rdf.org/pharmgkb_vocabulary:train_subset_final";
            List<Statement> trainStmts = new ArrayList<>();
            trainStmts = tripleStore.getStatementsFromStrings(null, trainURI, null);
            for (Statement stmt : trainStmts) {
                instances.add(stmt.getSubject());
                labels.add(stmt.getObject());
            }
            Map<Value, Double> labelMap = new HashMap<Value, Double>();
            double[] targets = new double[labels.size()];
            double t = 0;
            int m = 0;
            for (Value label : labels) {
                if (!labelMap.containsKey(label)) {
                    t += 1;
                    labelMap.put(label, t);
                }
                targets[m] = labelMap.get(label);
                m++;
            }
            Map<Double, Value> reslabelMap = EvaluationUtils.reverseLabelMap(labelMap);
            //double[] targets = createTargets(labels);
            String testURI = "http://bio2rdf.org/pharmgkb_vocabulary:test_set";
            List<Statement> testStmts = new ArrayList<>();
            testStmts = tripleStore.getStatementsFromStrings(null, testURI, null);
            for (Statement stmt : testStmts) {
                instances.add(stmt.getSubject());
                labels.add(stmt.getObject());
            }
            System.out.println("training instances: " + trainStmts.size());
            System.out.println("test instances: " + testStmts.size());
            RDFData data = new RDFData(tripleStore, instances, blackList);

            /// kernels
            String[] kernels = {"RDFGraphListWLSubTreeKernel", "RDFTreeWLSubTreeKernel", "RDFGraphListWalkCountKernel", "RDFTreeWalkCountKernel",
                "RDFGraphListWalkCountApproxKernelMkII", "RDFTreeWalkCountIDEQApproxKernelMkII", "RDFRootWalkCountKernel", "RDFRootWLSubTreeKernel",
                "RDFGraphListWLSubTreeApproxKernel", "RDFTreeWLSubTreeIDEQApproxKernel", "RDFRootWalkCountIDEQApproxKernel", "RDFRootWLSubTreeIDEQApproxKernel"};
             
            // kernels params
            // For details about the last 3 , see https://github.com/Data2Semantics/mustard/blob/master/mustard-kernels/src/main/java/org/data2semantics/mustard/kernels/graphkernels/singledtgraph/DTGraphWLSubTreeIDEQApproxKernel.java
            int[] pathLength_values = {0, 4, 6, 8, 10, 12, 15};
            int[] depth_values = {4, 6, 8, 10, 12, 15};
            boolean[] inference_values = {false};
            boolean[] normalize_values = {true};
            int[] minFreqs = {0,2,4,8,16,32};
            int[] maxPrevNBHs = {0, 1, 5};
            int[] maxLabelCards = {1, 10000};
            // svm params
            double[] cs = {1, 10, 100, 1000};
            //double[] nus = {0.1, 0.3, 0.5, 0.7, 0.9};
            LibSVMParameters svmParms = new LibSVMParameters(LibSVMParameters.C_SVC, cs);
            svmParms.setNumFolds(10);
            svmParms.setEvalFunction(new F1());
            svmParms.setProbEstimates(true);
            svmParms.setVerbosity(LibSVMParameters.VERBOSITY_DEFAULT);

            Map<Kernel, double[][]> currentKernels = new HashMap<>();
            for (String kernelName : kernels) {
                List<Object[]> paramsCombinations = createKernelParams(kernelName, pathLength_values, depth_values, inference_values, normalize_values, minFreqs, maxPrevNBHs, maxLabelCards);
                for (Object[] params : paramsCombinations) {
                    ComputeResult result = compute(kernelName, params, data);
                    currentKernels.put(result.getKernel(), result.getMatrix());
                }
            }
            //System.out.println("number of kernels: "+currentKernels.size());
            Map<Kernel, double[][]> trainKernels = new HashMap<>();
            Map<Kernel, double[][]> testKernels = new HashMap<>();
            for (Kernel k : currentKernels.keySet()) {
                //System.out.println("rows: "+currentKernels.get(k).length);
                //System.out.println("cols: "+currentKernels.get(k)[0].length);
                trainKernels.put(k, KernelUtils.trainSubset(currentKernels.get(k), 0, trainStmts.size()));
                testKernels.put(k, KernelUtils.testSubset(currentKernels.get(k), trainStmts.size(), trainStmts.size() + testStmts.size()));
            }
            LibSVMModel model = LibSVM.trainSVMModelWithMultipleKernels(trainKernels, targets, svmParms);
            long toc = System.currentTimeMillis();
            System.out.println("prediction kernel computaion time: " + (toc - tic));
            Prediction[] pred = LibSVM.testSVMModelWithMultipleKernels(model, testKernels);
            for (int e = 0; e < pred.length; e++) {
                Arrays.sort(pred, new Comparator<Prediction>() {
                    @Override
                    public int compare(Prediction p1, Prediction p2) {
                        return Double.valueOf(Math.max(p2.getDecisionValue()[0], p2.getDecisionValue()[1])).compareTo(Math.max(p1.getDecisionValue()[0], p1.getDecisionValue()[1]));
                    }
                });
                //Arrays.sort(pred);
                String instance = testStmts.get(pred[e].getIndex()).getSubject().stringValue();
                String label = reslabelMap.get(pred[e].getLabel()).stringValue();
                Double prob = Math.max(pred[e].getDecisionValue()[0], pred[e].getDecisionValue()[1]);
                System.out.println(instance + " : " + label + " : " + prob);
            }

        }

    }

    public static List<Object[]> createKernelParams(String kernelName, int[] pathLength_values, int[] depth_values, boolean[] inference_values, boolean[] normalize_values, int[] minFreqs, int[] maxPrevNBHs, int[] maxLabelCards) {
        List<Object[]> all = new ArrayList<>();
        switch (kernelName) {
            case "RDFGraphListWLSubTreeKernel":
                for (int depth : depth_values) {
                    for (int pathLength : pathLength_values) {
                        if (pathLength ==0 || pathLength == depth) {
                            for (boolean inference : inference_values) {
                                for (boolean normalize : normalize_values) {
                                    all.add(new Object[]{pathLength, depth, inference, true, true, normalize});
                                }
                            }
                        }
                    }
                }
                break;
            case "RDFTreeWLSubTreeKernel":
            case "RDFGraphListWalkCountKernel":
            case "RDFTreeWalkCountKernel":
                for (int depth : depth_values) {
                    for (int pathLength : pathLength_values) {
                        if (pathLength ==0 || pathLength == depth) {
                            for (boolean inference : inference_values) {
                                for (boolean normalize : normalize_values) {
                                    all.add(new Object[]{pathLength, depth, inference, normalize});
                                }
                            }
                        }
                    }
                }
                break;
            case "RDFGraphListWalkCountApproxKernelMkII":
            case "RDFTreeWalkCountIDEQApproxKernelMkII":
                for (int depth : depth_values) {
                    for (int pathLength : pathLength_values) {
                        if (pathLength ==0 || pathLength == depth) {
                            for (boolean inference : inference_values) {
                                for (int minFreq : minFreqs) {
                                    for (boolean normalize : normalize_values) {
                                        all.add(new Object[]{pathLength, depth, inference, minFreq, normalize});
                                    }
                                }
                            }
                        }
                    }
                }
                break;
            case "RDFRootWalkCountKernel":
            case "RDFRootWLSubTreeKernel":
                for (int pathLength : pathLength_values) {
                    for (boolean inference : inference_values) {
                        for (boolean normalize : normalize_values) {
                            all.add(new Object[]{pathLength, inference, normalize});
                        }
                    }
                }
                break;
            case "RDFGraphListWLSubTreeApproxKernel":
            case "RDFTreeWLSubTreeIDEQApproxKernel":
                for (int depth : depth_values) {
                    for (int pathLength : pathLength_values) {
                        if (pathLength ==0 || pathLength == depth) {
                            for (boolean inference : inference_values) {
                                for (boolean normalize : normalize_values) {
                                    all.add(new Object[]{pathLength, depth, inference, true, true, maxPrevNBHs, maxLabelCards, minFreqs, normalize});
                                }

                            }
                        }
                    }
                }
                break;
            case "RDFRootWalkCountIDEQApproxKernel":
                for (int pathLength : pathLength_values) {
                    for (boolean inference : inference_values) {
                        for (int minFreq : minFreqs) {
                            for (boolean normalize : normalize_values) {
                                all.add(new Object[]{pathLength, inference, minFreq, normalize});
                            }
                        }
                    }
                }
                break;
            case "RDFRootWLSubTreeIDEQApproxKernel":
                for (int pathLength : pathLength_values) {
                    for (boolean inference : inference_values) {
                        for (boolean normalize : normalize_values) {
                            all.add(new Object[]{pathLength, inference, maxPrevNBHs, maxLabelCards, minFreqs, normalize});
                        }
                    }
                }
                break;
            default:
                System.out.println("Invalid kernel name.");
                break;
        }
        return all;
    }

    public static ComputeResult compute(String kernelName, Object[] params, RDFData data) throws ClassNotFoundException, InstantiationException, IllegalAccessException, IllegalArgumentException, InvocationTargetException, NoSuchMethodException {

        Object kernel = Class.forName("org.data2semantics.mustard.kernels.graphkernels.rdfdata." + kernelName).getConstructors()[0].newInstance(params);
        Method computeFeatureVectors = kernel.getClass().getMethod("computeFeatureVectors", RDFData.class);
        Kernel k = (Kernel) kernel;
        SparseVector[] fvs = (SparseVector[]) computeFeatureVectors.invoke(kernel, data);
        double[][] matrix = KernelUtils.initMatrix(data.getInstances().size(), data.getInstances().size());
        matrix = KernelUtils.computeKernelMatrix(fvs, matrix);
        ComputeResult r = new ComputeResult(matrix, k);
        return r;
    }
}
