package hmc;

import beast.beagle.branchmodel.HomogeneousBranchModel;
import beast.beagle.sitemodel.GammaSiteRateModel;
import beast.beagle.treelikelihood.BeagleTreeLikelihood;
import beast.beagle.treelikelihood.PartialsRescalingScheme;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Patterns;
import beast.evolution.datatype.Nucleotides;
import beast.evolution.tree.Tree.MissingTaxonException;
import beast.evolution.util.Taxa;
import beast.evolution.util.Taxon;
import beast.evolution.util.Units.Type;
import beast.evomodel.branchratemodel.DefaultBranchRateModel;
import beast.evomodel.sitemodel.GammaSiteModel;
import beast.evomodel.sitemodel.SiteModel;
import beast.evomodel.speciation.BirthDeathModel;
import beast.evomodel.speciation.BirthDeathModel.TreeType;
import beast.evomodel.speciation.SpeciationLikelihood;
import beast.evomodel.substmodel.FrequencyModel;
import beast.evomodel.substmodel.HKY;
import beast.evomodel.tree.TreeModel;
import beast.inference.model.CompoundLikelihood;
import beast.inference.model.Likelihood;
import beast.inference.model.Parameter;
import beast.inference.model.Parameter.Default;
import beast.inference.operators.OperatorFailedException;
import beast.math.MathUtils;
import beast.xml.XMLObject;
import beast.xml.XMLParseException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Arrays;

/**
 * @author Arman Bilge
 */
public class Surface {

    private static final char[] alphabet = "abcdefghijklmnopqrstuvwxyz".toUpperCase().toCharArray();

    public static void main(final String... args) throws MissingTaxonException, OperatorFailedException, FileNotFoundException, XMLParseException {

        MathUtils.setSeed(666);

        final Taxa taxa = new Taxa();
        for (int i = 0; i < 3; ++i) {
            taxa.addTaxon(new Taxon("" + alphabet[i]));
        }
        final XMLObject xo = new XMLObject(null) {
            @Override
            public XMLObject getChild(String name) {
                if (name.equals("birthRate"))
                    return new XMLObject(null) {
                        @Override
                        public <T> T getChild(Class<T> c) {
                            return (T) new Parameter.Default(1.0);
                        }
                    };
                else
                    return new XMLObject(null) {
                        @Override
                        public <T> T getChild(Class<T> c) {
                            return (T) new Parameter.Default(0.0);
                        }
                    };
            }
            @Override
            public <T> T getChild(Class<T> c) {
                return (T) taxa;
            }
            @Override
            public boolean hasId() {
                return false;
            }
        };
        final TreeModel tree = new TreeModel(new BirthDeathSimulator().parseXMLObject(xo));
        final SiteModel site = new GammaSiteModel(new HKY(8.0, new FrequencyModel(Nucleotides.INSTANCE, new Default(new double[]{0.35, 0.30, 0.20, 0.15}))), new Default(0.1), null, 0, null);
        final Alignment alignment = new SequenceSimulator(tree, site, new DefaultBranchRateModel(), 512).simulate();
        final Likelihood like = new CompoundLikelihood(Arrays.<Likelihood>asList(
                new BeagleTreeLikelihood(new Patterns(alignment),
                                         tree,
                                         new HomogeneousBranchModel(new beast.beagle.substmodel.HKY(8, new beast.beagle.substmodel.FrequencyModel(Nucleotides.INSTANCE, new Default(new double[]{0.25, 0.25, 0.25, 0.25})))),
                                         new GammaSiteRateModel("", new Default(0.1), null, 0, null),
                                         new DefaultBranchRateModel(),
                                         null,
                                         false,
                                         PartialsRescalingScheme.NONE),
                new SpeciationLikelihood(tree, new BirthDeathModel(new Parameter.Default(1.0), null, null, null, TreeType.LABELED, Type.SUBSTITUTIONS), null)
        ));
        System.out.println(tree);
        final PrintWriter pw = new PrintWriter(new File("surface.dat"));
        pw.println("x y z");
        final int stepsize = 16;
        for (double x = 0/stepsize; x <= 2.0; x += 1.0/stepsize) {
            for (double y = 1.0/stepsize; y <= 2.0; y += 1.0/stepsize) {
                tree.setNodeHeight(tree.getInternalNode(1), x + y);
                tree.setNodeHeight(tree.getInternalNode(0), y);
                pw.println("" + x + " " + y + " " + -like.getLogLikelihood() + "");
            }
            pw.println();
        }
        pw.close();

//        final MyHamilton integrator = new MyHamilton(like, new Parameter[]{tree.createNodeHeightsParameter(true, true, false)}, null, Math.pow(2, -14), 1, 0.0, 1.0, CoercionMode.COERCION_OFF);
//        integrator.getP().adoptParameterValues(new Parameter.Default(new double[] {1.0, 1.0, 1.0}));
//        final JFrame frame = new JFrame();
//        final TreePane pane = new TreePane();
//        pane.setTreeLayout(new RectilinearTreeLayout());
//        System.out.println(pane.isTransformBranchesOn());
//        frame.add(pane);
//        frame.setBounds(0, 0, 500, 500);
//        frame.setVisible(true);
//        final PrintWriter pw = new PrintWriter(new File("data.txt"));
//        pw.println("x y z");
//        for (int i = 0; i < 3200000; ++i) {
//            if (i % 100 == 0) pw.println(String.join(" ", (Iterable<String>) IntStream.range(0, tree.getInternalNodeCount()).mapToObj(tree::getInternalNode).mapToDouble(tree::getNodeHeight).mapToObj(Double::toString)::iterator));
//            pane.setTree(tree.asJeblTree());
//            like.makeDirty();
//            integrator.doOperation();
//        }
//        pw.close();
    }


}
