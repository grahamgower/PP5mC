# Copyright (c) 2015,2017 Graham Gower <graham.gower@gmail.com>
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

import os, os.path
import jobscript
from graphviz import Digraph

class Pipeline:
    """
    For the creation and submission of a collection of related jobs.
    """

    class JobScript(jobscript.JobScript):

        def sub(self, afterok=None):
            """
            Submit the job to the queue, see jobscript.Jobscript.sub().

            Job IDs and the afterok dependencies are used to populate the
            graphviz dependency graph.

            A .success file is created when the job succeeds, which is
            used as a marker to avoid rerunning the job if the pipeline
            is restarted.
            This allows the entire pipeline to be restarted when it has
            failed part way through (and the problem fixed), maintaining
            dependencies yet avoiding rerunning of successful jobs.
            """

            filename = "{}/{}.sh".format(self.pipeline.acctdir, self.tmpl.name)
            success_filename = "{}.success".format(filename)

            if os.path.exists(success_filename):
                # if the success file already exists, do nothing
                return None

            # when the job succeeds, this file is created
            self.add_cmd("touch {}".format(success_filename))

            # write out the job script
            jobscript.JobScript.write(self, filename)

            if afterok:
                if not isinstance(afterok, list):
                    afterok = [afterok]
                else:
                    # filter Nones
                    afterok = [x for x in afterok if x is not None]

            if self.pipeline.fake_id > 0:
                job_id = str(self.pipeline.fake_id)
                self.pipeline.fake_id += 1
            else:
                job_id = jobscript.JobScript.sub(self, afterok)

            node_txt = "{}\nmem={}, cpus={}, time={}-{}:{}:{}".format(
                    self.tmpl.name, self.tmpl.mem, self.tmpl.cpus,
                    self.tmpl.days, self.tmpl.hours, self.tmpl.mins,
                    self.tmpl.secs)
            # update digraph
            self.pipeline.dot.node(job_id, node_txt)
            if afterok:
                for dep in afterok:
                    self.pipeline.dot.edge(dep, job_id)

            return job_id


    def __init__(self, queue, email, acctdir="./", dryrun=False):
        """
        Create a new pipeline on the specified queue.

        queue   -- passed to JobScript()
        email   -- passed to JobScript()
        acctdir -- folder for empty files that indicate job completion.
        dryrun -- if true, jobs will not be submitted. Fake job IDs will
                   be generated instead, for building the dependency graph.
        """
        self.queue = queue
        self.email = email
        self.acctdir = acctdir

        if dryrun:
            self.fake_id = 1
        else:
            self.fake_id = 0

        self.dot = Digraph()

    def new_job(self, name, res):
        """A wrapper for creating the JobScript."""
        j = self.JobScript(name=name, email=self.email, queue=self.queue, res=res)
        j.tmpl.output = "{}/{}.out".format(self.acctdir, name)
        j.pipeline = self
        return j

    def print_graph(self):
        """Output the dependency graph as produced by graphviz."""
        self.dot.render("dependency_graph")
