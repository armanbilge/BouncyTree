package hmc;

import beast.beagle.branchmodel.HomogeneousBranchModel;
import beast.beagle.sitemodel.GammaSiteRateModel;
import beast.beagle.treelikelihood.BeagleTreeLikelihood;
import beast.beagle.treelikelihood.PartialsRescalingScheme;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Patterns;
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
import java.util.Arrays;

/**
 * @author Arman Bilge
 */
public class BouncyTree {

    private static final char[] alphabet = "abcdefghijklmnopqrstuvwxyz".toUpperCase().toCharArray();

    public static void main(final String... args) throws MissingTaxonException, OperatorFailedException {

        final Taxa taxa = new Taxa();
        for (int i = 0; i < 16; ++i) {
            taxa.addTaxon(new Taxon("" + alphabet[i]));
        }
        final ConstantPopulation constant = new ConstantPopulation(Type.YEARS);
        constant.setN0(1);
        final TreeModel tree = new TreeModel(new CoalescentSimulator().simulateTree(taxa, constant));
        final SiteModel site = new GammaSiteModel(new HKY(8.0, new FrequencyModel(Nucleotides.INSTANCE, new Default(new double[]{0.25, 0.25, 0.25, 0.25}))));
        final Alignment alignment = new SequenceSimulator(tree, site, new DefaultBranchRateModel(), 2048).simulate();
        final Likelihood like = new CompoundLikelihood(Arrays.<Likelihood>asList(
                new BeagleTreeLikelihood(new Patterns(alignment),
                                         tree,
                                         new HomogeneousBranchModel(new beast.beagle.substmodel.HKY(8, new beast.beagle.substmodel.FrequencyModel(Nucleotides.INSTANCE, new Default(new double[]{0.25, 0.25, 0.25, 0.25})))),
                                         new GammaSiteRateModel(""),
                                         new DefaultBranchRateModel(),
                                         null,
                                         false,
                                         PartialsRescalingScheme.NONE),
                new CoalescentLikelihood(tree, null, null, new ConstantPopulationModel(new Default(1.0), Type.YEARS))
        ));

        final HamiltonUpdate integrator = new HamiltonUpdate(like, new Parameter[]{tree.createNodeHeightsParameter(false, true, false)}, null, Math.pow(2, -11), 1, 0.0, 1.0, CoercionMode.COERCION_OFF);

        final JFrame frame = new JFrame();
        final TreePane pane = new TreePane();
        pane.setTreeLayout(new RectilinearTreeLayout());
        frame.add(pane);
        frame.setBounds(0, 0, 500, 500);
        frame.setVisible(true);
        while (true) {
            pane.setTree(tree.asJeblTree());
            integrator.doOperation();
        }
    }

}
