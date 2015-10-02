/*
Copyright (c) 2015 Mengxiao Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
package com.sweetdum.hmm;

/**
 * HMMForwardBackwardManipulator is used to calculate alpha/beta/gamma array to
 * help the solving for HMM question.
 *
 * Created by Mengxiao Lin on 2015/10/2.
 */
public class HMMForwardBackwardManipulator {
    /**
     * the hmm representation bind to it.
     */
    private HMMRepresentation hmm;
    /**
     * the end time of the calculations
     */
    int endTime;
    /**
     * the observation sequence
     */
    int o[];
    /**
     * forward alpha array:
     * alpha[t][i] := P(o1 o2 ... ot-1, X[t]=s[i] | model)
     */
    double alpha[][];
    /**
     * backward beta array:
     * beta[t][i] := P (ot ot+1 ... oT | X[t]=s[i], model)
     */
    double beta[][];
    /**
     * gamma[t][i] := P (X[t]=s[i] | O,model)
     * is used to find the best state sequence
     */
    double gamma[][];

    /**
     * create a forward backward manipulator
     * @param hmm the HMMRepresentation bound to the manipulator
     * @param observations the observation ids array
     */
    public HMMForwardBackwardManipulator(HMMRepresentation hmm,int observations[]) {
        this.hmm = hmm;
        this.endTime = observations.length;
        this.o = observations;
        calculateAlpha();
        calculateBeta();
        calculateGamma();
    }

    /**
     * calculate alpha form the observation sequence
     */
    public void calculateAlpha() {
        //initial
        int N = hmm.stateCount;
        this.alpha = new double[endTime+1][N];
        for (int i = 0; i < N; ++i) {
            alpha[0][i] = hmm.pi[i];
        }
        //induction
        for (int t = 1; t <= endTime; ++t) {
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    alpha[t][i] += alpha[t - 1][j] * hmm.a[j][i] * hmm.b[j][i][o[t - 1]];
                }
            }
        }
    }
    /**
     * calculate beta form the observation sequence
     */
    public void calculateBeta() {
        //initial
        int N = hmm.stateCount;
        this.beta = new double[endTime + 1][N];
        for (int i = 0; i < N; ++i) beta[endTime][i] = 1;
        //induction
        for (int t = endTime - 1; t >= 0; --t) {
            for (int i = 0; i < N; ++i) {
                for (int j = 1; j < N; ++j) {
                    beta[t][i] += beta[t + 1][j] * hmm.a[i][j] * hmm.b[i][j][o[t]];
                }
            }
        }
    }
    /**
     * calculate gamma form alpha & beta
     * alpha & beta should have been prepared
     */
    public void calculateGamma(){
        /*
        use the formula gamma[t][i]= alpha[t][i]beta[t][i] /sigma(alpha[t][i],beta[t][i])
         */
        int N=hmm.stateCount;
        gamma=new double[endTime][N];
        for (int t=0;t<endTime;++t){
            double count=0;
            for (int i=0;i<N;++i){
                gamma[t][i]=alpha[t][i]*beta[t][i];
                count+=gamma[t][i];
            }
            for (int i=0;i<N;++i){
                gamma[t][i]/=count;
            }
        }
    }
}
