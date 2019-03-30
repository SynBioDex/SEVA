package SBOLpSEVA;


import java.io.File;
import java.net.URI;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.sbolstandard.core2.Component;
import org.sbolstandard.core2.ComponentDefinition;
import org.sbolstandard.core2.AccessType;
import org.sbolstandard.core2.RefinementType;
import org.sbolstandard.core2.SBOLReader;
import org.sbolstandard.core2.SBOLWriter;
import org.sbolstandard.core2.RestrictionType;
import org.sbolstandard.core2.SBOLDocument;
import org.sbolstandard.core2.OrientationType;
import org.sbolstandard.core2.Range;
import org.sbolstandard.core2.Sequence;
import org.sbolstandard.core2.SequenceAnnotation;
import org.sbolstandard.core2.SequenceConstraint;
import org.sbolstandard.core2.SequenceOntology;
import org.sbolstandard.core2.Collection;

public class Template3 extends SBOLDocument {
	public static void main(String[] args) throws Throwable {
//TODO
		String sevaURI="http://seva.cnb.csic.es/";
		
		String sevaPrefix="seva";	
		
		
		SBOLDocument document = new SBOLDocument();
		
		document.addNamespace(URI.create(sevaURI), sevaPrefix);
		document.setTypesInURIs(true);
		document.setDefaultURIprefix(sevaURI);
				
		
		Sequence seqPacI = document.createSequence(
				"seqPacI",
				"",
				"TTAATTAA", 
				URI.create(sevaURI)
				);
		
		Sequence seqSpeI = document.createSequence(
				"seqSpeI",
				"",
				"ACTAGT", 
				URI.create(sevaURI)
				);
		
		Sequence seqT0 = document.createSequence(
				"seqT0",
				"",
				"CTTGGACTCCTGTTGATAGATCCAGTAATGACCTCAGAACTCCATCTGGATT"
				+ "TGTTCAGAACGCTCGGTTGCCGCCGGGCGTTTTTTATTGGTGAGAATCCAG", 
				URI.create(sevaURI)
				);
		
		Sequence seqSanDI = document.createSequence(
				"seqSanDI",
				"",
				"GGGTCCC", 
				URI.create(sevaURI)
				);
		
		Sequence seqSwaI = document.createSequence(
				"seqSwaI",
				"",
				"ATTTAAAT", 
				URI.create(sevaURI)
				);
		
		Sequence seqPshAI = document.createSequence(
				"seqPshAI",
				"",
				"GACNNNNGTC", 
				URI.create(sevaURI)
				);
		
		Sequence seqOriT = document.createSequence(
				"seqOriT",
				"",
				"CTTTTCCGCTGCATAACCCTGCTTCGGGGTCATTATAGCGATTTTTTCG"
				+ "GTATATCCATCCTTTTTCGCACGATATACAGGATTTTGCCAAAGGGTTCGTGTAGA"
				+ "CTTTCCTTGGTGTATCCAACGGCGTCAGCCGGGCAGGATAGGTGAAGTAGGCCCAC"
				+ "CCGCGAGCGGGTGTTCCTTCTTCACTGTCCCTTATTCGCACCTGGCGGTGCTCAAC"
				+ "GGGAATCCTGCTCTGCGAGGCTGGCCGTA", 
				URI.create(sevaURI)
				);
		
		Sequence seqFseI = document.createSequence(
				"seqFseI",
				"",
				"GGCCGGCC", 
				URI.create(sevaURI)
				);
		
		Sequence seqAscI = document.createSequence(
				"seqAscI",
				"",
				"GGCGCGCC", 
				URI.create(sevaURI)
				);
		
		Sequence seqScarT0 = document.createSequence(
				"seqScarT0",
				"",
				"CAATAATTACG", 
				URI.create(sevaURI)
				);
		
		Sequence seqScarSmTc = document.createSequence(
				"seqScarSmTc",
				"",
				"ATTTACGT", 
				URI.create(sevaURI)
				);
		
		Sequence seqScarKmGm = document.createSequence(
				"seqScarKmGm",
				"",
				"CGCGCGTTGTC", 
				URI.create(sevaURI)
				);
		
		Sequence seqHok_Shok = document.createSequence(
				"seqHok_Shok",
				"",
				"AACAAACTCCGGGAGGCAGCGTGATGCGGCAACAATCACACGGAT"
				+ "TTCCCGTGAACGGTCTGAATGAGCGGATTATTTTCAGGGAAAGTGAGTGTGGTCAG"
				+ "CGTGCAGGTATATGGGCTATGATGTGCCCGGCGCTTGAGGCTTTCTGCCTCATGAC"
				+ "GTGAAGGTGGTTTGTTGCCGTGTTGTGTGGCAGAAAGAAGATAGCCCCGTAGTAAG"
				+ "TTAATTTTCATTAACCACCACGAGGCATCCCTATGTCTAGTCCACATCAGGATAGC"
				+ "CTCTTACCGCGCTTTGCGCAAGGAGAAGAAGGCCATGAAACTACCACGAAGTTCCC"
				+ "TTGTCTGGTGTGTGTTGATCGTGTGTCTCACACTGTTGATATTCACTTATCTGACA"
				+ "CGAAAATCGCTGTGCGAGATTCGTTACAGAGACGGACACAGGGAGGTGGCGGCTTT"
				+ "CATGGCTTACGAATCCGGTAAGTAGCAACCTGGAGGCGGGCGCAGGCCCGCCTTTT"
				+ "CAGGACTGATGCTGGTCTGACTACTGAAGCGCCTTTATAAAGGGGCTGCTGGTTCG"
				+ "CCGGTAGCCCCTTTCTCCTTGCTGATGTTGT", 
				URI.create(sevaURI)
				);
		
		Sequence seqT1 = document.createSequence(
				"seqT1",
				"",
				"CAGCTGTCTAGGGCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACA"
			+ "AACAACAGATAAAACGAAAGGCCCAGTCTTTCGACTGAGCCTTTCGTTTTATTTGAT"
			+ "GCCT", 
				URI.create(sevaURI)
				);
		
		
		ComponentDefinition pSeva = document.createComponentDefinition(
				"pSeva",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		
		
		ComponentDefinition PacI = document.createComponentDefinition(
				"PacI",
				"",	
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		
		PacI.setName("PacI");
		PacI.setDescription("Seva PacI");	
		PacI.addSequence(seqPacI.getIdentity());
		PacI.addRole(SequenceOntology.RESTRICTION_ENZYME_RECOGNITION_SITE);
		
		
		ComponentDefinition SpeI = document.createComponentDefinition(
				"SpeI",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		
		SpeI.setName("SpeI");
		SpeI.setDescription("Seva SpeI");	
		SpeI.addSequence(seqSpeI.getIdentity());
		SpeI.addRole(SequenceOntology.RESTRICTION_ENZYME_RECOGNITION_SITE);

		
		ComponentDefinition T0 = document.createComponentDefinition(
				"T0",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		
		T0.setName("T0");
		T0.setDescription("Seva T0");	
		T0.addSequence(seqT0.getIdentity());
		T0.addRole(SequenceOntology.TERMINATOR);

		
		ComponentDefinition SanDI = document.createComponentDefinition(
				"SanDI",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		
		SanDI.setName("SanDI");
		SanDI.setDescription("Seva SanDI");	
		SanDI.addSequence(seqSanDI.getIdentity());
		SanDI.addRole(SequenceOntology.RESTRICTION_ENZYME_RECOGNITION_SITE);

		
		ComponentDefinition SwaI = document.createComponentDefinition(
				"SwaI",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		
		SwaI.setName("SwaI");
		SwaI.setDescription("Seva SwaI");	
		SwaI.addSequence(seqSwaI.getIdentity());
		SwaI.addRole(SequenceOntology.RESTRICTION_ENZYME_RECOGNITION_SITE);

		
		ComponentDefinition PshAI = document.createComponentDefinition(
				"PshAI",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		
		PshAI.setName("PshAI");
		PshAI.setDescription("Seva PshAI");	
		PshAI.addSequence(seqPshAI.getIdentity());
		PshAI.addRole(SequenceOntology.RESTRICTION_ENZYME_RECOGNITION_SITE);

		
		ComponentDefinition OriT = document.createComponentDefinition(
				"OriT",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		
		OriT.setName("OriT");
		OriT.setDescription("Seva OriT");	
		OriT.addSequence(seqOriT.getIdentity());
		
		
		ComponentDefinition FseI = document.createComponentDefinition(
				"FseI",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
			
		FseI.setName("FseI");
		FseI.setDescription("Seva FseI");	;
		FseI.addSequence(seqFseI.getIdentity());
		FseI.addRole(SequenceOntology.RESTRICTION_ENZYME_RECOGNITION_SITE);

		
		ComponentDefinition AscI =document.createComponentDefinition(
				"AscI",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		
		AscI.setName("AscI");
		AscI.setDescription("Seva AscI");	
		AscI.addSequence(seqAscI.getIdentity());
		AscI.addRole(SequenceOntology.RESTRICTION_ENZYME_RECOGNITION_SITE);

		
		ComponentDefinition ScarT0 = document.createComponentDefinition(
				"ScarT0",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		
		ScarT0.setName("ScarT0");
		ScarT0.setDescription("Seva ScarT0");	
		ScarT0.addSequence(seqScarT0.getIdentity());
		ScarT0.addRole(SequenceOntology.TERMINATOR);

		
		ComponentDefinition ScarSmTc = document.createComponentDefinition(
				"ScarSmTc",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		
		ScarSmTc.setName("ScarSmTc");
		ScarSmTc.setDescription("Seva ScarSmTc");	
		ScarSmTc.addSequence(seqScarSmTc.getIdentity());
		ScarSmTc.addRole(SequenceOntology.GENE);

		
		ComponentDefinition ScarKmGm = document.createComponentDefinition(
				"ScarKmGm",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		
		ScarKmGm.setName("ScarKmGm");
		ScarKmGm.setDescription("ScarKmGm");	
		ScarKmGm.addSequence(seqScarKmGm.getIdentity());
		ScarKmGm.addRole(SequenceOntology.GENE);

		
		ComponentDefinition Hok_Shok = document.createComponentDefinition(
				"Hok_Shok",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		
		Hok_Shok.setName("Hok_Shok");
		Hok_Shok.setDescription("Hok_Shok");	
		Hok_Shok.addSequence(seqHok_Shok.getIdentity());
		Hok_Shok.addRole(SequenceOntology.ENGINEERED_GENE);

		
		ComponentDefinition Antibiotic = document.createComponentDefinition(
				"Antibiotic",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		Antibiotic.addRole(SequenceOntology.GENE);

		
		ComponentDefinition Cargo = document.createComponentDefinition(
				"Cargo",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		Cargo.addRole(SequenceOntology.ENGINEERED_GENE);

		
		ComponentDefinition OriV = document.createComponentDefinition(
				"OriV",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		OriV.addRole(SequenceOntology.ORIGIN_OF_REPLICATION);

		
		ComponentDefinition T1 = document.createComponentDefinition(
				"T1",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		
		T1.setName("T1");
		T1.setDescription("example T1");
		T1.addSequence(seqT1.getIdentity());
		T1.createSequenceAnnotation("ScarT1a", "ScarT1aloc", 1, 6, OrientationType.INLINE);
		T1.createSequenceAnnotation("T1seq", "T1seqloc", 7, 111, OrientationType.INLINE);
		T1.createSequenceAnnotation("ScarT1b", "ScarT1bloc", 112, 112, OrientationType.INLINE);
		T1.addRole(SequenceOntology.TERMINATOR);
		
		ComponentDefinition VariableRegion = document.createComponentDefinition(
				"VariableRegion",
				"",
				new HashSet<URI>(Arrays.asList(ComponentDefinition.DNA)));
		
		
		
		Component ComPacI = pSeva.createComponent("PacI", AccessType.PUBLIC, PacI.getIdentity());
    	//SequenceAnnotation anno1=pSeva.createSequenceAnnotation("anno1",  "loc1", start, end, OrientationType.INLINE);
		//anno1.setComponent(ComPacI.getIdentity());
		
		Component ComSpeI = pSeva.createComponent("SpeI", AccessType.PUBLIC, SpeI.getIdentity());

		Component ComT0 = pSeva.createComponent("T0", AccessType.PUBLIC, T0.getIdentity());
	
		Component ComSanDI = pSeva.createComponent("SanDI", AccessType.PUBLIC, SanDI.getIdentity());
	
		Component ComSwaI = pSeva.createComponent("SwaI", AccessType.PUBLIC, SwaI.getIdentity());
	
		Component ComPshAI = pSeva.createComponent("PshAI", AccessType.PUBLIC, PshAI.getIdentity());
	
		Component ComOriT = pSeva.createComponent("OriT", AccessType.PUBLIC, OriT.getIdentity());
		
		Component ComFseI = pSeva.createComponent("FseI", AccessType.PUBLIC, FseI.getIdentity());
	
		Component ComAscI = pSeva.createComponent("AscI", AccessType.PUBLIC, AscI.getIdentity());
	
		Component ComScarT0 = pSeva.createComponent("ScarT0", AccessType.PUBLIC, ScarT0.getIdentity());
		
		Component ComScarSmTc = pSeva.createComponent("ScarSmTc", AccessType.PUBLIC, ScarSmTc.getIdentity());
		
		Component ComScarKmGm = pSeva.createComponent("ScarKmGm", AccessType.PUBLIC, ScarKmGm.getIdentity());

		Component ComHok = pSeva.createComponent("Gadget", AccessType.PUBLIC, Hok_Shok.getIdentity());

		Component ComAntibiotic = pSeva.createComponent("Antibiotic", AccessType.PUBLIC, Antibiotic.getIdentity());
	
		Component ComCargo = pSeva.createComponent("Cargo", AccessType.PUBLIC, Cargo.getIdentity());
		
		Component ComOriV = pSeva.createComponent("OriV", AccessType.PUBLIC, OriV.getIdentity());

		Component ComT1 = pSeva.createComponent("T1", AccessType.PUBLIC, T1.getIdentity());
		
		Component ComVariableRegion = pSeva.createComponent("VariableRegion", AccessType.PUBLIC, VariableRegion.getIdentity());
		
		pSeva.createSequenceConstraint("PacICons", RestrictionType.PRECEDES, ComPacI.getIdentity(), ComCargo.getIdentity());
		pSeva.createSequenceConstraint("CargoCons", RestrictionType.PRECEDES, ComCargo.getIdentity(), ComSpeI.getIdentity());
		pSeva.createSequenceConstraint("SpeICons", RestrictionType.PRECEDES, ComSpeI.getIdentity(), ComT0.getIdentity());
		pSeva.createSequenceConstraint("T0Cons", RestrictionType.PRECEDES, ComT0.getIdentity(), ComSanDI.getIdentity());
		pSeva.createSequenceConstraint("SwaICons", RestrictionType.PRECEDES, ComSwaI.getIdentity(), ComAntibiotic.getIdentity());
		pSeva.createSequenceConstraint("OriTCons", RestrictionType.PRECEDES, ComOriT.getIdentity(), ComFseI.getIdentity());
		pSeva.createSequenceConstraint("FseICons", RestrictionType.PRECEDES, ComFseI.getIdentity(), ComOriV.getIdentity());
		pSeva.createSequenceConstraint("OriVCons", RestrictionType.PRECEDES, ComOriV.getIdentity(), ComAscI.getIdentity());
		pSeva.createSequenceConstraint("AscICons", RestrictionType.PRECEDES, ComAscI.getIdentity(), ComT1.getIdentity());
		pSeva.createSequenceConstraint("T1Cons", RestrictionType.PRECEDES, ComT1.getIdentity(), ComPacI.getIdentity());

		//TODO

		File file = new File ("/home/bioinfo/Escritorio/Template3.rdf");
		
		SBOLWriter.write(document,file);
		
	
	}
	
	
	
}