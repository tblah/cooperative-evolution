// reimplements https://eprints.soton.ac.uk/264277/
package CooperativeEvolution;

import scala.collection.mutable.ListBuffer;
import scala.util.Random;
import org.sameersingh.scalaplot._
import org.sameersingh.scalaplot.Implicits._
//import org.sameersingh.scalaplot.jfreegraph._
import org.sameersingh.scalaplot.gnuplot._

class Paper extends Common {
    // parameters from paper
    val growth_rate_cooperative = 0.018;
    val growth_rate_selfish = 0.02;
    val consumption_rate_cooperative = 0.1;
    val consumption_rate_selfish = 0.2;

    // stores statistics by generation for reproducing paper fig 2
    type PopType = Population;

    // draw graphs and save as right.png and left.png
    def draw_graphs: Unit = {
        val stats = previous_pops.toList;
        val x = breeze.linalg.linspace(0.0, number_of_generations, stats.length).toArray.toSeq;

        // left graph
        val left_data = new XYData();
        left_data += new MemXYSeries(x, stats.map(p => p.selfish / p.total).toSeq, "Selfish resource usage");
        left_data += new MemXYSeries(x, stats.map(p => p.big / p.total).toSeq, "Large group size");

        val left_chart = new XYChart("", left_data, x = Axis(label = "Generation"), y = Axis(label = "Global Frequency", range = Some((0.0, 0.8))));
        left_chart.showLegend = true;

        //val left_plotter = new JFGraphPlotter(left_chart);
        //left_plotter.gui();
        GnuplotPlotter.png(left_chart, "./", "left");

        // right graph
        val right_data = new XYData();
        right_data += new MemXYSeries(x, stats.map(p => p.cooperative_small / p.total).toSeq, "Cooperative + Small");
        right_data += new MemXYSeries(x, stats.map(p => p.cooperative_big / p.total).toSeq, "Cooperative + Large");
        right_data += new MemXYSeries(x, stats.map(p => p.selfish_small / p.total).toSeq, "Selfish + Small");
        right_data += new MemXYSeries(x, stats.map(p => p.selfish_big / p.total).toSeq, "Selfish + Large");

        val right_chart = new XYChart("", right_data, x = Axis(label = "Generation"), y = Axis(label = "Global Genotype Frquency", range = Some((0.0, 1.0))));
        right_chart.showLegend = true;

        //val right_plotter = new JFGraphPlotter(right_chart);
        //right_plotter.gui();
        GnuplotPlotter.png(right_chart, "./", "right");
    }

    // initialise the migrant pool with pop_size individuals
    protected def initialise: Population = {
        // all 4 genotypes start in equal proportion
        assert((pop_size % 4) == 0, "pop_size divisable by 4");
        val pop_each = pop_size/4;
        new Population(pop_each, pop_each, pop_each, pop_each)
    }

    // assign individuals from the migrant population to groups
    // returns each group as a population in an Iterable
    protected def assign_to_groups(migrant_pool: Population): Iterable[Population] = {
        var working_migrant_pool = migrant_pool; // scala doesn't let me have mutable function arguments :(
        var ret = ListBuffer.empty[Population];

        // small groups
        while ((working_migrant_pool.small >= n_small)) {
            var new_small_group = new Population(0, 0, 0, 0);
            // fill up the small group with random individuals
            for (_ <- 1 to n_small) {
                working_migrant_pool.random_small(new_small_group) match {
                    case (a, b) => {
                        new_small_group = a;
                        working_migrant_pool = b;
                    }
                }
            }
            //println("Made a new small group " + new_small_group.toString());

            ret += new_small_group;
        }

        // big groups
        while ((working_migrant_pool.big >= n_big)) {
            var new_big_group = new Population(0, 0, 0, 0);
            // fill up the large group with random individuals
            for (_ <- 1 to n_big) {
                working_migrant_pool.random_big(new_big_group) match {
                    case (a, b) => {
                        new_big_group = a;
                        working_migrant_pool = b;
                    }
                }
            }
            //println("Made a new big group " + new_big_group.toString());

            ret += new_big_group;
        }

        ret.asInstanceOf[Iterable[Population]]
    }

    protected def empty_pop = new Population(0, 0, 0, 0);

    // class to store representations of populations and groups as frequencies of individuals with each genotype
    protected class Population(val cooperative_small: Double, val cooperative_big: Double, val selfish_small: Double, val selfish_big: Double) extends AbstractPopulation {
        val rand = new Random(/*scala.compat.currentTime*/) // todo: random seed

        override def toString: String = "Population(" + cooperative_small + ", " + cooperative_big + ", " + selfish_small + ", " + selfish_big + ")";

        lazy val small = cooperative_small + selfish_small
        lazy val big = cooperative_big + selfish_big
        lazy val total = small + big
        lazy val selfish = selfish_small + selfish_big
        lazy val cooperative = cooperative_small + cooperative_big

        // removes one individual from the migrant pool and puts it into group
        // returns (group, migrant_pool)
        def random_small(group: Population): (Population, Population) = {
            val choice = small * rand.nextDouble();
            if (choice > cooperative_small) {
                (new Population(group.cooperative_small, group.cooperative_big, group.selfish_small + 1, group.selfish_big), new Population(cooperative_small, cooperative_big, selfish_small - 1, selfish_big))
            } else {
                (new Population(group.cooperative_small + 1, group.cooperative_big, group.selfish_small, group.selfish_big), new Population(cooperative_small - 1, cooperative_big, selfish_small, selfish_big))
            }
        }

        // works like random_small
        // (group, migrant_pool)
        def random_big(group: Population): (Population, Population) = {
            val choice = big * rand.nextDouble();
            if (choice > cooperative_big) {
                (new Population(group.cooperative_small, group.cooperative_big, group.selfish_small, group.selfish_big + 1), new Population(cooperative_small, cooperative_big, selfish_small, selfish_big - 1))
            } else {
                (new Population(group.cooperative_small, group.cooperative_big + 1, group.selfish_small, group.selfish_big), new Population(cooperative_small, cooperative_big - 1, selfish_small, selfish_big))
            }
        }

        def +(that: Population): Population = {
            new Population(cooperative_small + that.cooperative_small, cooperative_big + that.cooperative_big, selfish_small + that.selfish_small, selfish_big + that.selfish_big)
        }

        def reproduce: Population = {
            val denominator = (cooperative_small * consumption_rate_cooperative * growth_rate_cooperative) + (cooperative_big * consumption_rate_cooperative * growth_rate_cooperative) + (selfish_small * consumption_rate_selfish * growth_rate_selfish) + (selfish_big * consumption_rate_selfish * growth_rate_selfish);
            
            val new_cooperative_small = cooperative_small * (1 - death_rate + growth_rate_cooperative * r_small / denominator);
            val new_selfish_small = selfish_small * (1 - death_rate + growth_rate_selfish * r_small / denominator);
            val new_cooperative_big = cooperative_big * (1 - death_rate + growth_rate_cooperative * r_big / denominator);
            val new_selfish_big = selfish_big * (1 - death_rate + growth_rate_selfish * r_big / denominator);

            new Population(new_cooperative_small, new_cooperative_big, new_selfish_small, new_selfish_big);
        }

        // new Population rescaled back to pop_size individuals
        def rescaled: Population = {
            val factor = total / pop_size;

            new Population(cooperative_small/factor, cooperative_big/factor, selfish_small/factor, selfish_big/factor);
        }
    }
}

