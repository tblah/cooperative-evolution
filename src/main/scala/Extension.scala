// extends https://eprints.soton.ac.uk/264277/
package CooperativeEvolution;

import scala.util.Random;
import scala.collection.mutable.ArrayBuffer;

import org.sameersingh.scalaplot._
import org.sameersingh.scalaplot.Implicits._
import org.sameersingh.scalaplot.jfreegraph._
import org.sameersingh.scalaplot.gnuplot._

class Extension extends Common {
    // model parameters
    val std_dev = 0.5;
    val start_growth_rate = 19; 
    override val number_of_generations = 1000;
    override val n_small = 2;
    override val pop_size = 1000;
    override val r_big = 100;
    override val r_small = 125; 

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

    // moved out of population so we can override it
    // reproduce for one generation
    protected def p_reproduce(p: Population): Population = {
        val denominator = p.individuals.map(i => i.growth_rate * i.growth_rate).fold(0: Double)(_ + _);
        val relative_fitnesses = p.individuals.map(i => {
            val r = if (i.prefers_small) r_small else r_big;
            i.growth_rate * r / denominator
        });
        val fitness_scaling_factor = relative_fitnesses.fold(0: Double)(_ + _);
        val reproduction_probability = relative_fitnesses.map(_ / fitness_scaling_factor);
        val combined = p.individuals zip reproduction_probability;

        val new_pop = ArrayBuffer.empty[Individual];

        // fitness proportional selection
        for (_ <- 1 to combined.length) { // select combined.length individuals
            val choice = Random.nextDouble();
            var total: Double = 0;
            var chosen_one: Option[Individual] = None;

            for (i <- combined) { // for each potential parent
                if (chosen_one == None) { // if we haven't chosen a parent yet
                    total += i._2; 
                    if (total >= choice) {
                        chosen_one = Some(i._1);
                    }
                }
            }

            new_pop += chosen_one.get.mutate;
        }

        Population(new_pop.toIndexedSeq)
    }

    // class to store representations of populations and groups
    case class Population(val individuals: scala.collection.immutable.IndexedSeq[Individual]) extends AbstractPopulation {
        // "rescaled" is a name left over from the original paper. This is called once all gropus are returned to the migrant pool.
        // mutate the whole migrant pool with fitness proportional selection so that we get some competition between small and large groups
        def rescaled: Population = reproduce;

        // new population with only the small members
        lazy val smalls = Population(individuals.filter({case Individual(s, _) => s}));

        // new population with only the large members
        lazy val bigs = Population(individuals.filter({case Individual(s, _) => !s}));

        def reproduce = p_reproduce(this); 
        
        // union
        def +(other: Population): Population = Population(individuals ++ other.individuals)
    }

    // not used in paper
    // draws a greyscale image with each pixel as the growth rate of an individual. Y axis is individuals in populations, X axis is generations of populations
    def draw_image_plot(min: Double, max: Double, data: IndexedSeq[Population], title: String) = {
        val img_scale = breeze.plot.GradientPaintScale(min, max);

        val m = breeze.linalg.DenseMatrix.zeros[Double](data.head.individuals.length, data.length);
        for (i <- 0 to data.length - 1) {
            m(::, i) := breeze.linalg.DenseVector(data(i).individuals.map(_.growth_rate).toArray)
        }

        val fig = new breeze.plot.Figure(title, 1, 1);
        val sub = fig.subplot(0);
        sub.xlabel = "Generations";
        sub.ylabel = "Population";
        sub += breeze.plot.image(m, scale=img_scale);
    }

    // draws histograms figure
    def draw_hists(min: Double, max: Double, pop: Population, title: String) {
        val smalls_data = pop.smalls.individuals.map(_.growth_rate);
        val bigs_data = pop.bigs.individuals.map(_.growth_rate);

        val fig = new breeze.plot.Figure(title, 3, 1);

        if (smalls_data.length != 0) {
            val smalls = fig.subplot(1);
            smalls.title = "Small Groups";
            smalls.xlabel = "Growth Rate"
            smalls.ylabel = "Frequency"
            smalls += breeze.plot.hist(smalls_data, 100);
            smalls.xlim(min, max);
        }

        if (bigs_data.length != 0) {
            val bigs = fig.subplot(2);
            bigs.title = "Large Groups";
            bigs.xlabel = "Growth Rate"
            bigs.ylabel = "Frequency"
            bigs += breeze.plot.hist(pop.bigs.individuals.map(_.growth_rate), 100);
            bigs.xlim(min, max);
        } 

        val all = fig.subplot(0);
        all.title = "All";
        all.xlabel = "Growth Rate"
        all.ylabel = "Frequency"
        all += breeze.plot.hist(pop.individuals.map(_.growth_rate), 100);
        all.xlim(min, max);

        fig.refresh
        fig.saveas("./hists.png")
    }

    // draw graphs (unused in the paper)
    def draw_graphs = {
        val max_growth = previous_pops.map(_.individuals.map(_.growth_rate).max).max;
        val min_growth = previous_pops.map(_.individuals.map(_.growth_rate).min).min;
        println(min_growth + " <= growth_rate <= " + max_growth);

        // image plots
        draw_image_plot(min_growth, max_growth, previous_pops, "Growth Rate Locus Evolution");
        draw_image_plot(min_growth, max_growth, previous_pops.map(_.smalls), "Small");
        draw_image_plot(min_growth, max_growth, previous_pops.map(_.bigs), "Large");
        
        // histograms
        draw_hists(min_growth, max_growth, previous_pops.last, "Final Population");
    }

    // return a new empty population
    def empty_pop: Population = Population(scala.collection.immutable.IndexedSeq()); 

    // returns initial migrant pool
    // overriden in the version used in the paper
    def initialise: Population = {
        assert((pop_size % 2) == 0);
        val size = pop_size / 2;
        val dist = breeze.stats.distributions.Gaussian(start_growth_rate, std_dev);

        val new_pop = ArrayBuffer.empty[Individual];

        // small individuals
        for (_ <- 1 to size) {
            var sample = dist.sample;
            while (sample < 0) {
                sample = dist.sample;
            }
            new_pop += new Individual(true, sample);
        }

        // large individuals
        for (_ <- 1 to size) {
            var sample = dist.sample;
            while (sample < 0) {
                sample = dist.sample;
            }
            new_pop += new Individual(false, sample);
        }

        new Population(new_pop.toIndexedSeq)
    } 

    // assign individuals from the migrant population to groups
    // returns each group as a population in an Iterable
    def assign_to_groups(migrant_pool: Population): Iterable[Population] = {
        val groups = ArrayBuffer.empty[Population];

        // shuffles migrant_pool.smalls and then groups it into groups of size n_small
        val small_groups = Random.shuffle(migrant_pool.smalls.individuals).grouped(n_small);
        // don't remove individuals in incomplete groups from the population as this creates a bias against large groups
        groups ++= small_groups/*.filter(_.length == n_small)*/.map(Population(_));

        val big_groups = Random.shuffle(migrant_pool.bigs.individuals).grouped(n_big);
        groups ++= big_groups/*.filter(_.length == n_big)*/.map(Population(_));

        groups
    }
}

// overrides initialise() so that we start with cooperative + small, cooperative + large, selfish + small, selfigh + large
class MultimodalExtension extends Extension {
    // average initial growth rates
    val selfish_start = 10;
    val cooperative_start = 1;

    // returns initial migrant pool
    override def initialise: Population = {
        assert((pop_size % 4) == 0);
        val size = pop_size / 4;

        val new_pop = ArrayBuffer.empty[Individual];

        def make_individuals(small: Boolean, mean: Double) {
            val dist = breeze.stats.distributions.Gaussian(mean, std_dev);
            for (_ <- 1 to size) {
                var sample = dist.sample;
                while (sample < 0) {
                    sample = dist.sample;
                }
                new_pop += new Individual(small, sample);
            }
        }

        make_individuals(true, selfish_start);
        make_individuals(true, cooperative_start);
        make_individuals(false, selfish_start);
        make_individuals(false, cooperative_start);

        new Population(new_pop.toIndexedSeq)
    }
}