import com.sweetdum.hmm.HMMRepresentation;
import com.sweetdum.hmm.HMMSolver;

/**
 * Test program
 * Created by Mengxiao Lin on 2015/10/2.
 */
public class Tester {
    public static void main(String args[]){
        HMMRepresentation hmm=new HMMRepresentation(2,4);
        hmm.setA(0,0,0.1);
        hmm.setA(0,1,0.9);
        hmm.setA(1,0,0.1);
        hmm.setA(1,1,0.9);

        hmm.setB(0,0,0,1); hmm.setB(1,1,1,1);
        hmm.setB(0,1,2,1); hmm.setB(1,0,3,1);
        hmm.setPi(0,0.9);hmm.setPi(1,0.1);

        int observations[][]=new int[100][];
        int states[][]=new int[100][];
        //best stage test
        for (int i=0;i<100;++i){
            observations[i]=hmm.produceObservations(10);
            HMMSolver solver=new HMMSolver(hmm,observations[i]);
            states[i]=solver.getBestStateSequence();
        }
        //modify hmm
        HMMRepresentation hmm2=new HMMRepresentation(2,4);
        hmm2.uniformParameter();
        HMMRepresentation oldHmm;
        for (int iter=0;iter<1000;++iter){
            oldHmm=(HMMRepresentation)hmm2.clone();
            HMMSolver solver = new HMMSolver(hmm2, observations[0]);
            for (int i=0;i<100;++i){
                solver.readNewObservation(observations[i]);
                solver.estimateParameter();
            }
            double error=oldHmm.diffWithOtherHMM(hmm2);
            System.out.println("Iter "+iter+": ERROR="+error);
            if (error<1e-4) break;
        }
        hmm2.printParameter();
        for (int item=0;item<100;++item){
            HMMSolver solver=new HMMSolver(hmm2,observations[item]);
            int state[]=solver.getBestStateSequence();
            System.out.println("Data #"+item);
            System.out.println("Observation:");
            for (int i=0;i<observations[item].length;++i){
                System.out.print(observations[item][i]+" ");
            }
            System.out.println();
            System.out.println("State by model 1:");
            for (int i=0;i<states[item].length;++i){
                System.out.print(states[item][i]+" ");
            }
            System.out.println();
            System.out.println("State by model 2:");
            for (int i=0;i<state.length;++i){
                System.out.print(state[i]+" ");
            }
            System.out.println();
            System.out.println();
        }
    }
}
