__author__ = 'brian'

import glob
import os
import jinja2
import shutil
import sys
import errno

from utils import AnalyzedSample, Plot

"""
This maps the figures to explanations.  Not every figure may be used, but they are all explained here regardless

It maps the filenames to the descriptions.
"""
figure_description_mapping = {
    '1a.Indel_size_distribution_n_sequences.png':
        """ This figures shows the size distribution (in raw sequence counts) of the insertions
            and deletions in your amplicon(s).  The bar at zero represents
            sequences without INDELs, some of which may be unmodified; sequences
            with simple base substitutions preserving the length would also fall into that category.
        """,
    '1b.Indel_size_distribution_percentage.png':
        """ This figures shows the size distribution (in percentage of total reads) of the insertions
            and deletions in your amplicon(s).  The bar at zero represents
            sequences without INDELs, some of which may be unmodified; sequences
            with simple base substitutions preserving the length would also fall into that category.
        """,
    '2.Unmodified_NHEJ_HDR_pie_chart.png':
        """ This figure summarizes the distribution of reads that are unmodified, those that show incorporation
            of a donor sequence via the homology directed repair (HDR) pathway, and those that are likely repaired
            via nonhomoglous end joining (NHEJ).  Below the pie chart is a schematic representation of the location of the
            guide RNA (sgRNA) and the putative cut point of Cas9 on the amplicon. For example, this is approximately 3
            base pairs 5' of the PAM for SpCas9 (S. pyogenes Cas9)
        """,
    '2.Unmodified_NHEJ_pie_chart.png':
        """ This figure summarizes the distribution of reads that are unmodified and those that are likely repaired
            via nonhomoglous end joining (NHEJ); without a donor sequence, reliable HDR events cannot be deciphered.
            Below the pie chart is a schematic representation of the location of the
            guide RNA (sgRNA) and the putative cut point of Cas9 on the amplicon. For example, this is approximately 3
            base pairs 5' of the PAM for SpCas9 (S. pyogenes Cas9)
        """,
    '3.Insertion_Deletion_Substitutions_size_hist.png':
        """ This plot splits out the insertions and deletions into separate plots and also shows the distribution
        for the number of sequences with simple base substitutions.
        """,
    '4a.Combined_Insertion_Deletion_Substitution_Locations.png':
        """ This plot summarizes the location of any alterations present in the sequences; this includes insertions,
         deletions, and base substitutions.  The dotted line indicates the expected cut point(s) of the nuclease and the
         shaded regions show the location of the sgRNA within the amplicon sequence.  As stated by the CRISPResso authors,
         "only sequence positions directly adjacent to insertions or directly affected by deletions or substitutions are
         plotted."
        """,
    '4b.Insertion_Deletion_Substitution_Locations_NHEJ.png':
            """ This plot summarizes the location of any NHEJ alterations present in the sequences; this includes insertions,
         deletions, and base substitutions.  The dotted line indicates the expected cut point(s) of the nuclease and the
         shaded regions show the location of the sgRNA within the amplicon sequence.  As stated by the CRISPResso authors,
         "only sequence positions directly adjacent to insertions or directly affected by deletions or substitutions are
         plotted."
        """,
    '4c.Insertion_Deletion_Substitution_Locations_HDR.png':
        """ This plot summarizes the location of HDR alterations present in the sequences; this includes insertions,
         deletions, and base substitutions.  The dotted line indicates the expected cut point(s) of the nuclease and the
         shaded regions show the location of the sgRNA within the amplicon sequence.  As stated by the CRISPResso authors,
         "only sequence positions directly adjacent to insertions or directly affected by deletions or substitutions are
         plotted."
        """,
    '4d.Insertion_Deletion_Substitution_Locations_Mixed_HDR_NHEJ.png':
        """ This plot summarizes the location of mixed HDR+NHEJ alterations present in the sequences; this includes insertions,
         deletions, and base substitutions.  The dotted line indicates the expected cut point(s) of the nuclease and the
         shaded regions show the location of the sgRNA within the amplicon sequence.  As stated by the CRISPResso authors,
         "only sequence positions directly adjacent to insertions or directly affected by deletions or substitutions are
         plotted."
        """,
    '4e.Position_dependent_average_indel_size.png':
        """ This plot characterizes the size of insertions of deletions as they relate to their mapped position
        within the reference amplicon.""",
    '5.Frameshift_In-frame_mutations_pie_chart.png':
        """ If coding sequences are present (and supplied to the CRISPResso analysis software), the pie chart will
        show the distribution of frameshift events due to the alterations.  Below that is a schematic representation
          of the reference amplicon, the likely Cas9 cleavage position, and the location of the coding sequence.
        """,
    '6.Frameshift_In-frame_mutation_profiles.png':
        """This plot shows the distribution of both frameshift and in-frame events due to indels present in the
         sequences relative to the reference amplicon.
     """,
    '8.Potential_Splice_Sites_pie_chart.png':
        """ This shows the predicted impact on splice sites.  As stated by the authors, "potential splice sites modified
        refers to reads in which either of the two intronic positions adjacent to exon junctions are disrupted."  Thus,
        this does not consider the specific nature of the disruptions, merely their presence.
        """,
    '7.Insertion_Deletion_Substitution_Locations_Noncoding.png':
        """ This similarly shows the distribution of mutagenesis events, but shows only those exlusive to non-coding
        regions.  Again, as stated by the authors, "only sequence positions directly adjacent to insertions or directly affected by deletions or substitutions are
         plotted."
        """
}

CRISPRESSO_DIR_PREFIX = 'CRISPResso_on_'
LIBRARY_DIR = 'lib'


def main(project_dir):
    this_directory = os.path.dirname(os.path.realpath(__file__))

    report_dir = os.path.abspath(os.path.join(project_dir, 'html_report'))
    try:
        os.mkdir(report_dir)
    except OSError as ex:
        if ex.errno != errno.EEXIST:
            sys.exit(1)

    env = jinja2.Environment(loader=jinja2.FileSystemLoader(this_directory))
    template = env.get_template('template.html')
    path_pattern = os.path.join(project_dir, CRISPRESSO_DIR_PREFIX + '*')
    analysis_dirs = glob.glob(path_pattern)
    sample_list = []
    for i, d in enumerate(analysis_dirs):
        plot_paths = glob.glob(os.path.join(d, '*.png'))
        relative_plot_paths = [os.path.relpath(x, report_dir) for x in plot_paths]
        plot_list = [Plot(x, figure_description_mapping[os.path.basename(x)])
                     for x in relative_plot_paths]
        if i == 0:
            example_plots = plot_list
        sample_name = os.path.basename(d)[len(CRISPRESSO_DIR_PREFIX):]
        sample_list.append(AnalyzedSample(sample_name, plot_list))

    context = {'example_plots': example_plots, 'analyzed_samples': sample_list}

    completed_report_path = os.path.join(report_dir, 'report.html')
    with open(completed_report_path, 'w') as outfile:
        outfile.write(template.render(context))

    # copy the library files as well:
    libdir = os.path.join(this_directory, LIBRARY_DIR)
    destination = os.path.join(report_dir, LIBRARY_DIR)
    # shutil does not allow overwrite, so if the libraries were there,
    # need to rm them first
    if os.path.isdir(destination):
        shutil.rmtree(destination)
    shutil.copytree(libdir, destination)

project_dir = os.path.abspath(sys.argv[1])
main(project_dir)





