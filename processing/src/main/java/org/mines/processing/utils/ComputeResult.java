/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.mines.processing.utils;

import org.data2semantics.mustard.kernels.Kernel;

/**
 *
 * @author yassin
 */
public class ComputeResult {

        public ComputeResult(double[][] matrix, Kernel kernel) {
            this.matrix = matrix;
            this.kernel = kernel;
        }

        private double[][] matrix;

        public double[][] getMatrix() {
            return matrix;
        }

        private Kernel kernel;

        public Kernel getKernel() {
            return kernel;
        }

    }
