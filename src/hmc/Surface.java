package hmc;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Patterns;
import beast.evolution.datatype.Nucleotides;
import beast.evolution.tree.Tree.MissingTaxonException;
import beast.evolution.util.Taxa;
import beast.evolution.util.Taxon;
import beast.evomodel.branchratemodel.DefaultBranchRateModel;
import beast.evomodel.sitemodel.GammaSiteModel;
import beast.evomodel.sitemodel.SiteModel;
import beast.evomodel.substmodel.FrequencyModel;
import beast.evomodel.substmodel.HKY;
import beast.evomodel.tree.TreeModel;
import beast.evomodel.treelikelihood.TreeLikelihood;
import beast.inference.hamilton.HamiltonUpdate;
import beast.inference.hamilton.KineticEnergy.Fixed;
import beast.inference.loggers.ArrayLogFormatter;
import beast.inference.loggers.Logger;
import beast.inference.loggers.MCLogger;
import beast.inference.mcmc.MCMC;
import beast.inference.model.CompoundLikelihood;
import beast.inference.model.CompoundParameter;
import beast.inference.model.Likelihood;
import beast.inference.model.Parameter;
import beast.inference.model.Parameter.Default;
import beast.inference.operators.CoercionMode;
import beast.inference.operators.MCMCOperator;
import beast.inference.operators.OperatorFailedException;
import beast.inference.operators.ScaleOperator;
import beast.inference.operators.UniformOperator;
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
        final SiteModel site = new GammaSiteModel(new HKY(8.0, new FrequencyModel(Nucleotides.INSTANCE, new Default(new double[]{0.35, 0.30, 0.20, 0.15}))), new Default(1.0), null, 0, null);
        final Alignment alignment = new SequenceSimulator(tree, site, new DefaultBranchRateModel(), 512).simulate();
        final Likelihood like = new CompoundLikelihood(Arrays.<Likelihood>asList(
                new TreeLikelihood(new Patterns(alignment), tree, site, new DefaultBranchRateModel(), null, true, false, false, true, false)
//                new SpeciationLikelihood(tree, new BirthDeathModel(new Parameter.Default(1.0), null, null, null, TreeType.LABELED, Type.SUBSTITUTIONS), null)
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

        final MCMC mcmc = new MCMC("mcmc");

        tree.setNodeHeight(tree.getInternalNode(1), 2);
        tree.setNodeHeight(tree.getInternalNode(0), 1);

        {
            final ArrayLogFormatter formatter = new ArrayLogFormatter(true);
            final MCLogger logger = new MCLogger(formatter, 1, false);
            logger.add(like);
            logger.add(tree.createNodeHeightsParameter(true, true, false));
            mcmc.init(64, like, new MCMCOperator[]{new ScaleOperator(tree.getRootHeightParameter(), 0.5, CoercionMode.COERCION_OFF, 1.0), new UniformOperator(tree.createNodeHeightsParameter(false, true, false), 1)}, new Logger[]{logger});
//            mcmc.run();
        }

        tree.setNodeHeight(tree.getInternalNode(1), 2);
        tree.setNodeHeight(tree.getInternalNode(0), 1);

        {
            final ArrayLogFormatter formatter = new ArrayLogFormatter(true);
            final MCLogger logger = new MCLogger(formatter, 1, false);
            logger.add(like);
            logger.add(tree.createNodeHeightsParameter(true, true, false));
            mcmc.init(64, like, new MCMCOperator[]{new HamiltonUpdate(like, (CompoundParameter) tree.createNodeHeightsParameter(true, true, false), new Fixed(new double[][]{{1.0, 0}, {0, 1.0}}), 1.0/128, 4, 0.5, 1.0, CoercionMode.COERCION_OFF)}, new Logger[]{logger});
            mcmc.run();
        }

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
