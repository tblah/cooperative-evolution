package CooperativeEvolution

import scala.collection._

import org.sameersingh.scalaplot._
import org.sameersingh.scalaplot.Implicits._
import org.sameersingh.scalaplot.jfreegraph._
import org.sameersingh.scalaplot.gnuplot._

object Main extends App {
    def extension_avg = {
        val num_runs = 20;
        val e = new MultimodalExtension;

        def xy(data: immutable.IndexedSeq[immutable.IndexedSeq[Double]], title: String): MemXYSeries = {
            println("Averaging")
            val avgs = data.map((v: immutable.IndexedSeq[Double]) => v.fold(0.0)(_ + _) / v.length);
            val x = breeze.linalg.linspace(0.0, e.number_of_generations, avgs.length).toArray.toSeq;

            println("MemXYSeries")
            new MemXYSeries(x, avgs.toSeq, title);
        }

        def growth_graph(data: immutable.IndexedSeq[e.Population], title: String) = {
            val xydata = new XYData();
            xydata += xy(data.map(_.individuals.map(_.growth_rate)), "All");
            xydata += xy(data.map(_.bigs.individuals.map(_.growth_rate)), "Large Groups");
            xydata += xy(data.map(_.smalls.individuals.map(_.growth_rate)), "Small Groups");
            
            println("Plotting");
            val chart = new XYChart(title, xydata, x = Axis(label = "Generation"), y = Axis(label = "Average Growth Rate"));
            chart.showLegend = true;

            //(new JFGraphPlotter(chart)).gui();
            GnuplotPlotter.png(chart, "./", "growth-rates");
        }

        def draw_prop_large(data: IndexedSeq[IndexedSeq[e.Population]]) = {
            val x = breeze.linalg.linspace(0.0, e.number_of_generations, data.head.length).toArray.toSeq;
            val prop_large = new XYData();
    
            for (i <- data) {
                prop_large += new MemXYSeries(x, i.map(p => p.bigs.individuals.length.toDouble / p.individuals.length.toDouble).toSeq);
            }
    
            val prop_large_chart = new XYChart("", prop_large, x = Axis(label = "Generation"), y = Axis(label = "Proportion of individuals prefering a large group"));
            //(new JFGraphPlotter(prop_large_chart)).gui();
            GnuplotPlotter.png(prop_large_chart, "./", "prop-large");
        }

        val results = mutable.ArrayBuffer.empty[immutable.IndexedSeq[e.Population]];

        // do first run separately
        e.itterate;
        results += e.previous_pops.toIndexedSeq;
        //e.draw_graphs;
        //growth_graph(results.head, "First run");

        println("Calculating more results...");

        // get more results
        for (i <- 2 to num_runs) {
            println(i + "/" + num_runs);
            e.itterate
            results += e.previous_pops.toIndexedSeq
        }

        draw_prop_large(results);

        // average results
        println("joining populations...")
        val avgs = results.transpose.map(l => l.fold(e.empty_pop)(_ + _)).toIndexedSeq;

        growth_graph(avgs, "Average accross " + num_runs + " runs");

        println("Histograms");
        val min = avgs.last.individuals.map(_.growth_rate).min
        val max = avgs.last.individuals.map(_.growth_rate).max
        e.draw_hists(min, max, avgs.last, "Final Population Average accross " + num_runs + " runs");
        println("Done");
    }

    
    (new Paper).run;
    extension_avg
}
