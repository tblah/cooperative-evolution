// reimplements https://eprints.soton.ac.uk/264277/
package CooperativeEvolution;

import scala.collection.AbstractIterable;
import scala.collection.mutable.ListBuffer;
import scala.util.Random;
import breeze.plot._

class Paper {
    // parameters from paper
    val generations_in_group = 4; // t
    val growth_rate_cooperative = 0.018;
    val growth_rate_selfish = 0.02;
    val consumption_rate_cooperative = 0.1;
    val consumption_rate_selfish = 0.2;
    val pop_size = 4000;
    val number_of_generations = 120 * generations_in_group; // the graphs in the paper seem to only show data from when groups are dispersed into the migrant pool. The graphs go up to 120
    val death_rate = 0.1;
    val n_small = 4;
    val n_big = 40;
    val r_small = 4;
    val r_big = 50;

    // stores statistics by generation for reproducing paper fig 2
    private val previous_pops = ListBuffer.empty[Population];

    // run the simulation and display statistics
    def run: Unit = {
        assert((number_of_generations % generations_in_group) == 0, "N divisible by t");
        var groups = assign_to_groups(initialise);

        for (_ <- 1 to (number_of_generations / generations_in_group)) {
            val migrant_pool: Population = reproduce_within_groups(groups).rescaled;
            groups = assign_to_groups(migrant_pool);
        }

        val stats = previous_pops.toList;
        println("stats.length = " + stats.length);

        val figure = new Figure("My Implementation of Figure 2", 1, 2);
        val left = figure.subplot(0);
        left.xlabel = "Generation"
        left.ylabel = "Global Frequency"
        left.legend = true
        left.ylim(0, 0.8);
        val right = figure.subplot(1);
        right.xlabel = "Generation"
        right.ylabel = "Global Genotype Frequency"
        right.legend = true
        right.ylim(0, 1);

        val x = breeze.linalg.linspace(0.0, number_of_generations, stats.length);

        left += plot(x, stats.map(p => p.big / p.total), style = '.', name = "Large Group Size");
        left += plot(x, stats.map(p => p.selfish / p.total), name = "Selfish Resource Usage");

        right += plot(x, stats.map(p => p.cooperative_small / p.total), name = "Cooperative + Small");
        right += plot(x, stats.map(p => p.cooperative_big / p.total), style = '.', name = "Cooperative + Large");
        right += plot(x, stats.map(p => p.selfish_small / p.total), name = "Selfish + Small");
        right += plot(x, stats.map(p => p.selfish_big / p.total), style = '.', name = "Selfish + Large");

        figure.saveas("fig2.png");
    }

    // initialise the migrant pool with pop_size individuals
    private def initialise(): Population = {
        // all 4 genotypes start in equal proportion
        assert((pop_size % 4) == 0, "pop_size divisable by 4");
        val pop_each = pop_size/4;
        new Population(pop_each, pop_each, pop_each, pop_each)
    }

    // assign individuals from the migrant population to groups
    private def assign_to_groups(migrant_pool: Population): Iterable[Population] = {
        var working_migrant_pool = migrant_pool;
        var ret = ListBuffer.empty[Population];

        // small groups
        while ((working_migrant_pool.small >= n_small)) {
            var new_small_group = new Population(0, 0, 0, 0);
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

    // reproduce within groups for generations_in_group timesteps
    private def reproduce_within_groups(groups: Iterable[Population]): Population = {
        var working_groups = groups;

        for (_ <- 1 to generations_in_group) {
            // record statistics
            previous_pops += working_groups.fold(new Population(0, 0, 0, 0))(_ + _);

            // reproduce
            working_groups = working_groups.map(p => p.reproduce);
        }

        working_groups.fold(new Population(0, 0, 0, 0))(_ + _);
    }

    // class to store representations of populations and subpopulations as frequencies of individuals with each genotype
    private class Population(val cooperative_small: Double, val cooperative_big: Double, val selfish_small: Double, val selfish_big: Double) {
        val rand = new Random(/*scala.compat.currentTime*/)

        override def toString: String = "Paper.Population(" + cooperative_small.toString + ", " + cooperative_big.toString + ", " + selfish_small.toString + ", " + selfish_big.toString + ")";

        lazy val small = cooperative_small + selfish_small
        lazy val big = cooperative_big + selfish_big
        lazy val total = small + big
        lazy val selfish = selfish_small + selfish_big
        lazy val cooperative = cooperative_small + cooperative_big

        // (new, migrant_pool)
        def random_small(group: Population): (Population, Population) = {
            val choice = small * rand.nextDouble();
            if (choice > cooperative_small) {
                (new Population(group.cooperative_small, group.cooperative_big, group.selfish_small + 1, group.selfish_big), new Population(cooperative_small, cooperative_big, selfish_small - 1, selfish_big))
            } else {
                (new Population(group.cooperative_small + 1, group.cooperative_big, group.selfish_small, group.selfish_big), new Population(cooperative_small - 1, cooperative_big, selfish_small, selfish_big))
            }
        }

        // (new, migrant_pool)
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
            
            val new_cooperative_small = cooperative_small * (1 - death_rate + ((growth_rate_cooperative * r_small) / denominator));
            val new_selfish_small = selfish_small * (1 - death_rate + ((growth_rate_selfish * r_small) / denominator));
            val new_cooperative_big = cooperative_big * (1 - death_rate + ((growth_rate_cooperative * r_big) / denominator));
            val new_selfish_big = selfish_big * (1 - death_rate + ((growth_rate_selfish * r_big) / denominator));

            new Population(new_cooperative_small, new_cooperative_big, new_selfish_small, new_selfish_big);
        }

        // new Population rescaled back to N individuals
        def rescaled: Population = {
            val sum = total;
            val factor = sum / pop_size;

            val ret = new Population(cooperative_small/factor, cooperative_big/factor, selfish_small/factor, selfish_big/factor);

            //println("Population rescaled from " + sum + " to " + ret.total + ": " + ret.toString);
            ret
        }
    }
}

