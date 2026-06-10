/* 
 * Copyright 2016 Lutz Fischer <l.fischer@ed.ac.uk>.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package rappsilber.ms.crosslinker;

import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import rappsilber.config.ConfigurationParserException;
import rappsilber.config.RunConfig;
import rappsilber.ms.sequence.AminoAcid;
import rappsilber.ms.sequence.Peptide;
import rappsilber.ms.sequence.ions.loss.CleavableCrossLinkerPeptide;

/**
 * base class for crosslinker, that link to a set of given amino acids
 * @author Lutz Fischer <l.fischer@ed.ac.uk>
 */
public abstract class AminoAcidRestrictedCrossLinker extends CrossLinker{

    protected static class StubTarget {
        private final AminoAcid aminoAcid;
        private final boolean nTerminal;
        private final boolean cTerminal;

        protected StubTarget(AminoAcid aminoAcid, boolean nTerminal, boolean cTerminal) {
            this.aminoAcid = aminoAcid;
            this.nTerminal = nTerminal;
            this.cTerminal = cTerminal;
        }

        protected boolean matches(Peptide pep, int position) {
            if (aminoAcid != null && pep.nonLabeledAminoAcidAt(position) == aminoAcid) {
                return true;
            }
            if (nTerminal && position == 0 && pep.isNTerminal()) {
                return true;
            }
            return cTerminal && position == pep.length() - 1 && pep.isCTerminal();
        }
    }

    protected static class StubCondition {
        private final int sourceSite;
        private final StubTarget target;

        protected StubCondition(int sourceSite, StubTarget target) {
            this.sourceSite = sourceSite;
            this.target = target;
        }

        protected boolean matches(Peptide moiety0Peptide, int moiety0Position, Peptide moiety1Peptide, int moiety1Position) {
            if (sourceSite == 0) {
                return target.matches(moiety0Peptide, moiety0Position);
            }
            return moiety1Peptide != null && target.matches(moiety1Peptide, moiety1Position);
        }
    }

    /** all amino-acids, that can be cross-linked */
    protected HashMap<AminoAcid,Double> m_linkable;
    /** can the cross-linker react with the N-terminal */
    protected boolean            m_NTerminal = true;
    protected double             m_NTerminalWeight = 1;
    /** can the cross-linker react with the C-terminal */
    protected boolean            m_CTerminal = false;
    protected double             m_CTerminalWeight = 1;
    protected boolean            m_linksEverything = false;
    protected double             m_defaultAminoAcidWeight = Double.POSITIVE_INFINITY;
    protected HashSet<String> m_globalStubs = new HashSet<String>();
    protected HashMap<Integer, HashSet<String>> m_unconditionalStubs = new HashMap<Integer, HashSet<String>>();
    protected HashMap<Integer, HashMap<String, ArrayList<StubCondition>>> m_conditionalStubs = new HashMap<Integer, HashMap<String, ArrayList<StubCondition>>>();

    /**
     * creates and registers a new crosslinker (- definition)
     * @param Name
     * @param BaseMass
     * @param CrossLinkedMass
     * @param linkableAminoAcids
     */
    public AminoAcidRestrictedCrossLinker (String Name, double BaseMass, double CrossLinkedMass, HashSet<AminoAcid> linkableAminoAcids) {
        super(Name, BaseMass, CrossLinkedMass);
        m_linkable = new HashMap<AminoAcid, Double>();
        for (AminoAcid aa : linkableAminoAcids) {
            m_linkable.put(aa, 1.0);
            
            // I use this as part of the workaround for nary cross-links
            //m_linkable.put(AminoAcid.X, 1.0);
        }
    }    

    /**
     * creates and registers a new crosslinker (- definition)
     * @param Name
     * @param BaseMass
     * @param CrossLinkedMass
     * @param linkableAminoAcids
     */
    public AminoAcidRestrictedCrossLinker (String Name, double BaseMass, double CrossLinkedMass, HashMap<AminoAcid,Double> linkableAminoAcids) {
        super(Name, BaseMass, CrossLinkedMass);
        m_linkable = linkableAminoAcids;
    }

    /**
     * creates and registers a new crosslinker (- definition)
     * @param Name
     * @param BaseMass
     * @param CrossLinkedMass
     * @param linkableAminoAcids
     */
    public AminoAcidRestrictedCrossLinker (String Name, double BaseMass, double CrossLinkedMass, AminoAcid[] linkableAminoAcids) {
        super(Name, BaseMass, CrossLinkedMass);
        HashMap<AminoAcid,Double> linkable = new HashMap<AminoAcid,Double>(linkableAminoAcids.length);

        for (AminoAcid aa : linkableAminoAcids) {
            linkable.put(aa, 1.0);
        }

        m_linkable = linkable;
    }


    public double getAminoAcidWeight(AminoAcid AA) {
        Double w = 0d;
        if (!m_linkable.isEmpty()) {
            w = m_linkable.get(AA);
        } else if (m_linksEverything) {
            w = m_defaultAminoAcidWeight;
        }      
        if (w == null && m_linksEverything) {
            return m_defaultAminoAcidWeight;
        }
        return w == null ? Double.POSITIVE_INFINITY : w;
    }


    public void setlinksNTerm(boolean linksNTerm) {
        m_NTerminal = linksNTerm;
    }
    
    public void setNTermWeight(double w) {
        m_NTerminalWeight = w;
    }

    public void setlinksCTerm(boolean linksCTerm) {
        m_CTerminal = linksCTerm;
    }

    public void setCTermWeight(double w) {
        m_CTerminalWeight = w;
    }

    @Override
    public double getWeight(Peptide pep, int position) {

        double aaw = getAminoAcidWeight(pep.nonLabeledAminoAcidAt(position));

        if (position==0 && m_NTerminal && pep.isNTerminal()) {
            return Math.min(m_NTerminalWeight, aaw);
        }

        if (position==pep.length()-1 && m_CTerminal && pep.isCTerminal()) {
            return Math.min(m_CTerminalWeight, aaw);
        }

        return aaw;

    }

    @Override
    public double getWeight(Peptide pep, int position, int site) {
        return getWeight(pep, position);
    }

    @Override
    public boolean canProduceStub(Peptide stubPeptide, int stubSite, Peptide otherPeptide, int otherSite, String stubName) {
        return hasStubTarget(0, stubPeptide, stubSite, otherPeptide, otherSite, stubName);
    }

    @Override
    public boolean canProduceStub(String stubName) {
        if (m_globalStubs.contains(stubName)) {
            return true;
        }
        for (HashSet<String> siteStubs : m_unconditionalStubs.values()) {
            if (siteStubs.contains(stubName)) {
                return true;
            }
        }
        for (HashMap<String, ArrayList<StubCondition>> siteStubs : m_conditionalStubs.values()) {
            if (siteStubs.containsKey(stubName)) {
                return true;
            }
        }
        return false;
    }
    
    
    public Set<AminoAcid> getAASpecificity(int site) {
        return m_linkable.keySet();
    }

    @Override
    public boolean linksCTerminal(int site) {
        return m_CTerminal;
    }

    @Override
    public boolean linksNTerminal(int site) {
        return m_NTerminal;
    }

    protected void setLinksEverything(boolean linksEverything, double defaultWeight) {
        m_linksEverything = linksEverything;
        m_defaultAminoAcidWeight = defaultWeight;
    }

    protected void registerStub(String stubName) {
        m_globalStubs.add(stubName);
    }

    protected void registerStub(int site, String stubName) {
        HashSet<String> siteStubs = m_unconditionalStubs.get(site);
        if (siteStubs == null) {
            siteStubs = new HashSet<String>();
            m_unconditionalStubs.put(site, siteStubs);
        }
        siteStubs.add(stubName);
    }

    protected void registerConditionalStub(int site, int sourceSite, String stubName, StubTarget target) {
        HashMap<String, ArrayList<StubCondition>> siteStubs = m_conditionalStubs.get(site);
        if (siteStubs == null) {
            siteStubs = new HashMap<String, ArrayList<StubCondition>>();
            m_conditionalStubs.put(site, siteStubs);
        }
        ArrayList<StubCondition> targets = siteStubs.get(stubName);
        if (targets == null) {
            targets = new ArrayList<StubCondition>();
            siteStubs.put(stubName, targets);
        }
        targets.add(new StubCondition(sourceSite, target));
    }

    protected boolean hasStubTarget(int site, Peptide pep, int position, Peptide otherPeptide, int otherPosition, String stubName) {
        return hasStubTargetForMoieties(site, pep, position, otherPeptide, otherPosition, stubName);
    }

    protected boolean hasStubTargetForMoieties(int site, Peptide moiety0Peptide, int moiety0Position, Peptide moiety1Peptide, int moiety1Position, String stubName) {
        if (m_globalStubs.contains(stubName)) {
            return true;
        }
        HashSet<String> siteUnconditionalStubs = m_unconditionalStubs.get(site);
        if (siteUnconditionalStubs != null && siteUnconditionalStubs.contains(stubName)) {
            return true;
        }
        HashMap<String, ArrayList<StubCondition>> siteStubs = m_conditionalStubs.get(site);
        if (siteStubs == null) {
            return false;
        }
        ArrayList<StubCondition> targets = siteStubs.get(stubName);
        if (targets == null) {
            return false;
        }
        for (StubCondition target : targets) {
            if (target.matches(moiety0Peptide, moiety0Position, moiety1Peptide, moiety1Position)) {
                return true;
            }
        }
        return false;
    }

    protected static StubTarget parseStubTarget(String targetDefinition, RunConfig config) throws ConfigurationParserException {
        String target = targetDefinition.trim();
        String normalized = target.toLowerCase().replaceAll("-", "");
        if (normalized.contentEquals("nterm")) {
            return new StubTarget(null, true, false);
        }
        if (normalized.contentEquals("cterm")) {
            return new StubTarget(null, false, true);
        }
        AminoAcid aa = config.getAminoAcid(target);
        if (aa == null) {
            throw new ConfigurationParserException("Unknown conditional stub target: " + targetDefinition);
        }
        return new StubTarget(aa, false, false);
    }

    protected static void registerStubProducer(String stubName, double stubMass, RunConfig config) throws ParseException {
        CleavableCrossLinkerPeptide.parseArgs("MASS:" + stubMass + ";NAME:" + stubName, config);
    }

    protected static int parseStubConditionSourceSite(String targetDefinition, int defaultSite) {
        String normalized = targetDefinition.toUpperCase();
        if (normalized.startsWith("FIRSTLINKEDAMINOACIDS(") || normalized.startsWith("FIRST(")) {
            return 0;
        }
        if (normalized.startsWith("SECONDLINKEDAMINOACIDS(") || normalized.startsWith("SECOND(")) {
            return 1;
        }
        return defaultSite;
    }

    protected static String parseStubConditionTargets(String targetDefinition) {
        String target = targetDefinition.trim();
        int open = target.indexOf('(');
        if (open >= 0 && target.endsWith(")")) {
            return target.substring(open + 1, target.length() - 1);
        }
        return target;
    }

    protected static void parseConditionalStubList(String definition, int site, AminoAcidRestrictedCrossLinker crosslinker, RunConfig config) throws ConfigurationParserException, ParseException {
        for (String mapping : definition.split("\\|")) {
            String entry = mapping.trim();
            if (entry.isEmpty()) {
                continue;
            }
            String[] parts = entry.split("=", 2);
            if (parts.length != 2) {
                throw new ConfigurationParserException("Conditional stubs need the format target=stubName,mass,... but got " + entry);
            }
            int sourceSite = parseStubConditionSourceSite(parts[0].trim(), site);
            String[] targetDefinitions = parseStubConditionTargets(parts[0]).split(",");
            String[] stubs = parts[1].split(",");
            if (stubs.length % 2 != 0) {
                throw new ConfigurationParserException("Conditional stubs must be pairs of name and mass in " + entry);
            }
            for (int s = 0; s < stubs.length; s += 2) {
                String stubName = stubs[s].trim();
                double stubMass = Double.parseDouble(stubs[s + 1].trim());
                registerStubProducer(stubName, stubMass, config);
                for (String targetDefinition : targetDefinitions) {
                    StubTarget target = parseStubTarget(targetDefinition, config);
                    crosslinker.registerConditionalStub(site, sourceSite, stubName, target);
                }
            }
        }
    }

}
