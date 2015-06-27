package hmc;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Patterns;
import beast.evolution.datatype.Nucleotides;
import beast.evolution.tree.Tree.MissingTaxonException;
import beast.evolution.util.Taxa;
import beast.evolution.util.Taxon;
import beast.evolution.util.Units;
import beast.evomodel.branchratemodel.DefaultBranchRateModel;
import beast.evomodel.sitemodel.GammaSiteModel;
import beast.evomodel.sitemodel.SiteModel;
import beast.evomodel.speciation.BirthDeathModel;
import beast.evomodel.speciation.SpeciationLikelihood;
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
import beast.inference.operators.*;
import beast.inference.trace.Trace;
import beast.math.MathUtils;
import beast.xml.XMLObject;
import beast.xml.XMLParseException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;

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
                new TreeLikelihood(new Patterns(alignment), tree, site, new DefaultBranchRateModel(), null, true, false, false, true, false),
                new SpeciationLikelihood(tree, new BirthDeathModel(new Parameter.Default(1.0), null, null, null, BirthDeathModel.TreeType.LABELED, Units.Type.SUBSTITUTIONS), null)
        ));
        System.out.println(tree);

        {
            final PrintWriter pw = new PrintWriter(new File("surface.dat"));
            pw.println("x y z");
            final int stepsize = 32;
            for (double x = 8.0/stepsize; x <= 1.5; x += 1.0/stepsize) {
                for (double y = 1.0/stepsize; y <= 1.125; y += 1.0/stepsize) {
                    tree.setNodeHeight(tree.getInternalNode(1), x + y);
                    tree.setNodeHeight(tree.getInternalNode(0), y);
                    pw.println("" + x + " " + y + " " + -like.getLogLikelihood() + "");
                }
                pw.println();
            }
            pw.close();
        }

        final MCMC mcmc = new MCMC("mcmc");

        tree.setNodeHeight(tree.getInternalNode(1), 1.5);
        tree.setNodeHeight(tree.getInternalNode(0), 1);

        {
            final ArrayLogFormatter formatter = new ArrayLogFormatter(true);
            final MCLogger logger = new MCLogger(formatter, 1, false);
            logger.add(like);
            logger.add(tree.getRootHeightParameter());
            logger.add(tree.createNodeHeightsParameter(false, true, false));
            mcmc.init(64, like, new MCMCOperator[]{new ScaleOperator(tree.getRootHeightParameter(), 0.5, CoercionMode.COERCION_OFF, 1.0), new UniformOperator(tree.createNodeHeightsParameter(false, true, false), 1)}, new Logger[]{logger});
            mcmc.run();
            final PrintWriter pw = new PrintWriter(new File("mcmc.dat"));
            pw.println("x y z");
            List<Trace> traces = formatter.getTraces();
            for (int i = 0; i < traces.get(0).getValuesSize(); ++i) {
                final double y = (double) traces.get(3).getValue(i);
                final double x = (double) traces.get(2).getValue(i) - y;
                final double z = -(double) traces.get(1).getValue(i);
                pw.println(x + " " + y + " " + z);
            }
            pw.close();
        }

        tree.setNodeHeight(tree.getInternalNode(1), 1.5);
        tree.setNodeHeight(tree.getInternalNode(0), 1);

        {
            final ArrayLogFormatter formatter = new ArrayLogFormatter(true);
            final MCLogger logger = new MCLogger(formatter, 1, false);
            logger.add(like);
            logger.add(tree.getRootHeightParameter());
            logger.add(tree.createNodeHeightsParameter(false, true, false));
            mcmc.init(16, like, new MCMCOperator[]{new HamiltonUpdate(like, (CompoundParameter) tree.createNodeHeightsParameter(true, true, false), new Fixed(new double[][]{{1.0, 0}, {0, 1.0}}), 1.0/32, 4, 0.5, 1.0, CoercionMode.COERCION_OFF)}, new Logger[]{logger});
            mcmc.run();
            final PrintWriter pw = new PrintWriter(new File("hmc.dat"));
            pw.println("x y z");
            List<Trace> traces = formatter.getTraces();
            for (int i = 0; i < traces.get(0).getValuesSize(); ++i) {
                final double y = (double) traces.get(3).getValue(i);
                final double x = (double) traces.get(2).getValue(i) - y;
                final double z = -(double) traces.get(1).getValue(i);
                pw.println(x + " " + y + " " + z);
            }
            pw.close();

        }

    }


}
