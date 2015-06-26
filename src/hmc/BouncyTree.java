package hmc;

import beast.evolution.alignment.Alignment;
import beast.evolution.coalescent.CoalescentSimulator;
import beast.evolution.coalescent.ConstantPopulation;
import beast.evolution.datatype.Nucleotides;
import beast.evolution.tree.Tree.MissingTaxonException;
import beast.evolution.util.Taxa;
import beast.evolution.util.Taxon;
import beast.evolution.util.Units.Type;
import beast.evomodel.branchratemodel.DefaultBranchRateModel;
import beast.evomodel.coalescent.CoalescentLikelihood;
import beast.evomodel.coalescent.ConstantPopulationModel;
import beast.evomodel.sitemodel.GammaSiteModel;
import beast.evomodel.sitemodel.SiteModel;
import beast.evomodel.substmodel.FrequencyModel;
import beast.evomodel.substmodel.HKY;
import beast.evomodel.tree.TreeModel;
import beast.inference.hamilton.HamiltonUpdate;
import beast.inference.model.CompoundLikelihood;
import beast.inference.model.Likelihood;
import beast.inference.model.Parameter;
import beast.inference.model.Parameter.Default;
import beast.inference.operators.CoercionMode;
import beast.inference.operators.OperatorFailedException;
import figtree.treeviewer.TreePane;
import figtree.treeviewer.treelayouts.RectilinearTreeLayout;

import javax.swing.JFrame;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * @author Arman Bilge
 */
public class BouncyTree {

    private static final char[] alphabet = "abcdefghijklmnopqrstuvwxyz".toUpperCase().toCharArray();

    public static class MyHamilton extends HamiltonUpdate {
        public MyHamilton(Likelihood U, Parameter[] parameters, double[] massAttribute, double epsilon, int L, double alpha, @OperatorWeightAttribute double weight, CoercionMode mode) {
            super(U, parameters, massAttribute, epsilon, L, alpha, weight, mode);
        }

        public Parameter getP() {
            return p;
        }
    }

    public static void main(final String... args) throws MissingTaxonException, OperatorFailedException, FileNotFoundException {

        final Taxa taxa = new Taxa();
        for (int i = 0; i < 4; ++i) {
            taxa.addTaxon(new Taxon("" + alphabet[i]));
        }
        final ConstantPopulation constant = new ConstantPopulation(Type.YEARS);
        constant.setN0(6);
        final TreeModel tree = new TreeModel(new CoalescentSimulator().simulateTree(taxa, constant));
        tree.setNodeHeight(tree.getInternalNode(0), 1.0);
        tree.setNodeHeight(tree.getInternalNode(1), 3.0);
        tree.setNodeHeight(tree.getInternalNode(2), 9.0);
        final SiteModel site = new GammaSiteModel(new HKY(8.0, new FrequencyModel(Nucleotides.INSTANCE, new Default(new double[]{0.25, 0.25, 0.25, 0.25}))));
        final Alignment alignment = new SequenceSimulator(tree, site, new DefaultBranchRateModel(), 32).simulate();
        final Likelihood like = new CompoundLikelihood(Arrays.<Likelihood>asList(
//                new BeagleTreeLikelihood(new Patterns(alignment),
//                                         tree,
//                                         new HomogeneousBranchModel(new beast.beagle.substmodel.HKY(8, new beast.beagle.substmodel.FrequencyModel(Nucleotides.INSTANCE, new Default(new double[]{0.25, 0.25, 0.25, 0.25})))),
//                                         new GammaSiteRateModel(""),
//                                         new DefaultBranchRateModel(),
//                                         null,
//                                         false,
//                                         PartialsRescalingScheme.NONE),
                new CoalescentLikelihood(tree, null, null, new ConstantPopulationModel(new Default(6.0), Type.YEARS))
        ));

        final MyHamilton integrator = new MyHamilton(like, new Parameter[]{tree.createNodeHeightsParameter(true, true, false)}, null, Math.pow(2, -14), 1, 0.0, 1.0, CoercionMode.COERCION_OFF);
        integrator.getP().adoptParameterValues(new Parameter.Default(new double[] {1.0, 1.0, 1.0}));
        final JFrame frame = new JFrame();
        final TreePane pane = new TreePane();
        pane.setTreeLayout(new RectilinearTreeLayout());
        System.out.println(pane.isTransformBranchesOn());
        frame.add(pane);
        frame.setBounds(0, 0, 500, 500);
        frame.setVisible(true);
        final PrintWriter pw = new PrintWriter(new File("data.txt"));
        pw.println("x y z");
        for (int i = 0; i < 3200000; ++i) {
            if (i % 100 == 0) pw.println(String.join(" ", (Iterable<String>) IntStream.range(0, tree.getInternalNodeCount()).mapToObj(tree::getInternalNode).mapToDouble(tree::getNodeHeight).mapToObj(Double::toString)::iterator));
            pane.setTree(tree.asJeblTree());
            like.makeDirty();
            integrator.doOperation();
        }
        pw.close();
    }

}
