// reimplements https://eprints.soton.ac.uk/264277/
package CooperativeEvolution;

import scala.collection.AbstractIterable;
import scala.collection.mutable.ListBuffer;
import scala.util.Random;

class Paper {
    // parameters from paper
    val growth_rate_cooperative = 0.018;
    val growth_rate_selfish = 0.02;
    val consumption_rate_cooperative = 0.1;
    val consumption_rate_selfish = 0.2;
    val pop_size = 4000;
    val number_of_generations = 1000;
    val death_rate = 0.1;
    val n_small = 4;
    val n_big = 40;
    val r_small = 10;
    val r_big = 50;
    val generations_in_group = 4;

    // stores statistics by generation for reproducing paper fig 2
    private val previous_pops = ListBuffer.empty[Population];

    // run the simulation
    def run: Unit = {
        assert((number_of_generations % generations_in_group) == 0, "N divisible by t");
        var groups = assign_to_groups(initialise);

        for (_ <- 1 to (number_of_generations / generations_in_group)) {
            val migrant_pool: Population = reproduce_within_groups(groups).rescaled;
            groups = assign_to_groups(migrant_pool);
        }
    }

    // initialise the migrant pool with pop_size individuals
    private def initialise(): Population = {
        // all 4 genotypes start in equal proportion
        assert((pop_size % 4) == 0, "pop_size divisable by 4");
        val pop_each = pop_size/4;
        new Population(pop_each, pop_each, pop_each, pop_each)
    }

    // assign individuals from the migrant population to groups
    private def assign_to_groups(migrant_pool: Population): AbstractIterable[Population] = {
        var ret = ListBuffer.empty[Population];

        // small groups
        while ((migrant_pool.small >= n_small)) {
            val new_small_group = new Population(0, 0, 0, 0);
            for (_ <- 1 to n_small) {
                if (migrant_pool.random_is_selfish_small)
                    new_small_group.selfish_small += 1;
                else
                    new_small_group.cooperative_small += 1;
            }
            //println("Made a new small group " + new_small_group.toString());

            ret += new_small_group;
        }

        // big groups
        while ((migrant_pool.big >= n_big)) {
            val new_big_group = new Population(0, 0, 0, 0);
            for (_ <- 1 to n_big) {
                if (migrant_pool.random_is_selfish_big)
                    new_big_group.selfish_big += 1;
                else
                    new_big_group.cooperative_big += 1;
            }
            //println("Made a new big group " + new_big_group.toString());

            ret += new_big_group;
        }

        ret.asInstanceOf[AbstractIterable[Population]]
    }

    // reproduce within groups for generations_in_group timesteps
    private def reproduce_within_groups(groups: AbstractIterable[Population]): Population = {
        for (_ <- 1 to generations_in_group) {
            // record statistics
            previous_pops += groups.fold(new Population(0, 0, 0, 0))(_ + _);

            // reproduce
            groups.foreach(p => p.reproduce);
        }

        groups.fold(new Population(0, 0, 0, 0))(_ + _);
    }

    // class to store representations of populations and subpopulations as frequencies of individuals with each genotype
    private class Population(var cooperative_small: Double, var cooperative_big: Double, var selfish_small: Double, var selfish_big: Double) {
        val rand = new Random(/*scala.compat.currentTime*/)

        override def toString: String = "Paper.Population(" + cooperative_small.toString + ", " + cooperative_big.toString + ", " + selfish_small.toString + ", " + selfish_big.toString + ")";

        def small = cooperative_small + selfish_small
        def big = cooperative_big + selfish_big
        def total = small + big

        def random_is_selfish_small: Boolean = {
            val choice = small * rand.nextDouble();
            if (choice > cooperative_small) {
                selfish_small -= 1;
                true
            } else {
                cooperative_small -= 1;
                false
            }
        }

        def random_is_selfish_big: Boolean = {
            val choice = big * rand.nextDouble();
            if (choice > cooperative_big) {
                selfish_big -= 1;
                true
            } else {
                cooperative_big -= 1;
                false
            }
        }

        def +(that: Population): Population = {
            new Population(cooperative_small + that.cooperative_small, cooperative_big + that.cooperative_big, selfish_small + that.selfish_small, selfish_big + that.selfish_big)
        }

        def reproduce = {
            val denominator = cooperative_small * consumption_rate_cooperative * growth_rate_cooperative + cooperative_big * consumption_rate_cooperative * growth_rate_cooperative + selfish_small * consumption_rate_selfish * growth_rate_selfish + selfish_big * consumption_rate_selfish * growth_rate_selfish;
            
            cooperative_small = cooperative_small * (1 - death_rate + growth_rate_cooperative * r_small / denominator);
            selfish_small = selfish_small * (1 - death_rate + growth_rate_selfish * r_small / denominator);
            cooperative_big = cooperative_big * (1 - death_rate + growth_rate_cooperative * r_big / denominator);
            selfish_big = selfish_big * (1 - death_rate + growth_rate_selfish * r_big / denominator);
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

