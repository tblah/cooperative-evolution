// extends https://eprints.soton.ac.uk/264277/
package CooperativeEvolution;

import scala.collection.mutable.ListBuffer;
import scala.util.Random;

import org.sameersingh.scalaplot._
import org.sameersingh.scalaplot.Implicits._
import org.sameersingh.scalaplot.jfreegraph._
import org.sameersingh.scalaplot.gnuplot._

class Extension extends Common {
    // model parameters
    val consumption_scaling_factor: Double = 10;
    val std_dev = 0.005;
    val start_growth_rate = 0.019; 
    //override val number_of_generations = 1000;

    protected type PopType = Population;

    // individuals within the population
    case class Individual(val prefers_small: Boolean, val growth_rate: Double) {
        assert(growth_rate >= 0);
        def mutate = {
            val gaussian = breeze.stats.distributions.Gaussian(growth_rate, std_dev);

            // don't grow negatively
            var sample = gaussian.sample;
            while (sample < 0) {
                sample = gaussian.sample;
            }

            new Individual(prefers_small, sample)
        }
    }

    // class to store representations of populations and groups
    case class Population(val individuals: List[Individual]) extends AbstractPopulation {
        // not applicable to this population model because rescaling is done during reproduction
        def rescaled: Population = this;

        // population with only the small members
        lazy val smalls = Population(individuals.filter({case Individual(s, _) => s}));

        // population with only the large members
        lazy val bigs = Population(individuals.filter({case Individual(s, _) => !s}));

        // reproduce for one generation
        def reproduce: Population = {
            val denominator = individuals.map(i => i.growth_rate * i.growth_rate * consumption_scaling_factor).fold(0: Double)(_ + _);
            val relative_fitnesses = individuals.map(i => {
                val r = if (i.prefers_small) r_small else r_big;
                /*1 - death_rate +*/ i.growth_rate * r / denominator
            });
            val fitness_scaling_factor = relative_fitnesses.fold(0: Double)(_ + _);
            val reproduction_probability = relative_fitnesses.map(_ / fitness_scaling_factor);
            val combined = individuals zip reproduction_probability;

            val new_pop = ListBuffer.empty[Individual];

            // fitness proportional selection
            for (_ <- 1 to combined.length) { // select combined.length individuals
                val choice = Random.nextDouble();
                var total: Double = 0;
                var chosen_one: Option[Individual] = None;

                for (i <- combined) {
                    if (chosen_one == None)
                    total += i._2;               
                    if (total >= choice) {
                        chosen_one = Some(i._1);
                    }
                }

                new_pop += chosen_one.get.mutate;
            }

            Population(new_pop.toList)
        }
        
        // union
        def +(other: Population): Population = Population(individuals ++ other.individuals)
    }

    // draw graphs and save as right.png and left.png
    def draw_graphs = {
        val stats = previous_pops.toList;
//        val x = breeze.linalg.linspace(0.0, number_of_generations, stats.length).toArray.toSeq;

        val max_growth = stats.last.individuals.map(_.growth_rate).max;
        val min_growth = stats.last.individuals.map(_.growth_rate).min;
        println(min_growth + " <= growth_rate <= " + max_growth);

        // proportion large
        /*val prop_large = new XYData();
        prop_large += new MemXYSeries(x, stats.map(p => p.bigs.individuals.length.toDouble / p.individuals.length.toDouble).toSeq, "Large group size");

        val prop_large_chart = new XYChart("Proportion with large group allele", prop_large, x = Axis(label = "Generation"), y = Axis(label = "Frequency"));
        (new JFGraphPlotter(prop_large_chart)).gui();*/

        // image plots
        val img_scale = breeze.plot.GradientPaintScale(min_growth, max_growth);
        val m = breeze.linalg.DenseMatrix.zeros[Double](pop_size, stats.length); 
        for (i <- 0 to stats.length - 1) {
            m(::, i) := breeze.linalg.DenseVector(stats(i).individuals.map(_.growth_rate).toArray)
        }

        val fig2 = new breeze.plot.Figure("All - Growth rate by generation", 1, 1);
        val sub = fig2.subplot(0);
        sub.xlabel = "Generations";
        sub.ylabel = "Growth rate allele"
        sub += breeze.plot.image(m, scale=img_scale);

        fig2.saveas("extention-all.png");
        fig2.refresh

        val m2 = breeze.linalg.DenseMatrix.zeros[Double](pop_size/2, stats.length); 
        for (i <- 0 to stats.length - 1) {
            m2(::, i) := breeze.linalg.DenseVector(stats(i).smalls.individuals.map(_.growth_rate).toArray)
        }

        val fig21 = new breeze.plot.Figure("Small - Growth rate by generation", 1, 1);
        val sub1 = fig21.subplot(0);
        sub1.xlabel = "Generations";
        sub1.ylabel = "Growth rate allele"
        sub1 += breeze.plot.image(m2, scale=img_scale);

        fig21.saveas("extention-small.png");
        fig21.refresh

        for (i <- 0 to stats.length - 1) {
            m2(::, i) := breeze.linalg.DenseVector(stats(i).bigs.individuals.map(_.growth_rate).toArray)
        }

        val fig22 = new breeze.plot.Figure("Large - Growth rate by generation", 1, 1);
        val sub2 = fig22.subplot(0);
        sub2.xlabel = "Generations";
        sub2.ylabel = "Growth rate allele"
        sub2 += breeze.plot.image(m2, scale=img_scale);

        fig22.saveas("extention-big.png");
        fig22.refresh

        // histograms
        val fig3 = new breeze.plot.Figure("Final population histograms", 3, 1);

        val bigs = fig3.subplot(0);
        bigs.title = "Large groups histogram";
        bigs += breeze.plot.hist(stats.last.bigs.individuals.map(_.growth_rate), 100);
        bigs.xlim(min_growth, max_growth);

        val smalls = fig3.subplot(1);
        smalls.title = "Small groups histogram";
        smalls += breeze.plot.hist(stats.last.smalls.individuals.map(_.growth_rate), 100);
        smalls.xlim(min_growth, max_growth);

        val all = fig3.subplot(2);
        all.title = "All histogram";
        all += breeze.plot.hist(stats.last.individuals.map(_.growth_rate), 100);
        all.xlim(min_growth, max_growth);

        fig3.refresh
    }

    // return a new empty population
    def empty_pop: Population = Population(List()); 

    // returns initial migrant pool
    def initialise: Population = {
        assert((pop_size % 2) == 0);
        val size = pop_size / 2;
        val dist = breeze.stats.distributions.Gaussian(start_growth_rate, std_dev);

        val new_pop = ListBuffer.empty[Individual];

        // small individuals
        for (_ <- 1 to size) {
            /*var sample = dist.sample;
            while (sample < 0) {
                sample = dist.sample;
            }
            new_pop += new Individual(true, sample);*/
            new_pop += new Individual(true, start_growth_rate);
        }

        // large individuals
        for (_ <- 1 to size) {
            /*var sample = dist.sample;
            while (sample < 0) {
                sample = dist.sample;
            }
            new_pop += new Individual(false, sample);*/
            new_pop += new Individual(false, start_growth_rate);
        }

        new Population(new_pop.toList)
    } 

    // assign individuals from the migrant population to groups
    // returns each group as a population in an Iterable
    def assign_to_groups(migrant_pool: Population): Iterable[Population] = {
        val groups = ListBuffer.empty[Population];

        // shuffles migrant_pool.smalls and then groups it into groups of size n_small
        val small_groups = Random.shuffle(migrant_pool.smalls.individuals).grouped(n_small);
        // remove incomplete groups and wrap in Population() beofore adding to list of groups
        groups ++= small_groups.filter(_.length == n_small).map(Population(_));

        val big_groups = Random.shuffle(migrant_pool.bigs.individuals).grouped(n_big);
        groups ++= big_groups.filter(_.length == n_big).map(Population(_));

        groups
    }
}

