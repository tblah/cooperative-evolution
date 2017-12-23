// reimplements https://eprints.soton.ac.uk/264277/
package CooperativeEvolution;

import scala.collection.mutable.ListBuffer;

abstract class Common {
    // parameters from paper
    val generations_in_group = 4; // t
    val pop_size = 4000;
    val number_of_generations = 120 * generations_in_group; // the graphs in the paper seem to only show data from when groups are dispersed into the migrant pool. The graphs go up to 120
    val death_rate = 0.1;
    val n_small = 4;
    val n_big = 40;
    val r_small = 4;
    val r_big = 50;

    // stores statistics by generation for reproducing paper fig 2
    protected type PopType <: AbstractPopulation;
    val previous_pops = ListBuffer.empty[PopType];

    // class to store representations of populations and groups
    protected abstract class AbstractPopulation() {
        // scaled back to pop_size
        def rescaled: PopType; 

        // reproduce for one generation
        def reproduce: PopType;
        
        // union
        def +(other: PopType): PopType;
    }

    // itterate through each generation
    def itterate: Unit = {
        assert((number_of_generations % generations_in_group) == 0, "N divisible by t");

        // empty stats data if it is not already empty
        previous_pops.trimStart(previous_pops.length);

        // create the initial migrant pool and assign groups
        var groups: Iterable[PopType] = assign_to_groups(initialise);

        // for each group lifetime
        for (_ <- 1 to (number_of_generations / generations_in_group)) {
            // reproduce within groups and rescale population
            val migrant_pool = reproduce_within_groups(groups).rescaled;
            // assign the new migrant pool into new groups
            groups = assign_to_groups(migrant_pool);
        }
    }

    // reproduce within groups for generations_in_group timesteps
    private def reproduce_within_groups(groups: Iterable[PopType]): PopType = {
        var working_groups = groups; // can't have mutable function arguments

        for (_ <- 1 to generations_in_group) {
            // record statistics as a single population consisting of the sumation of all the groups
            previous_pops += working_groups.fold(empty_pop)(_ + _);

            // reproduce
            working_groups = working_groups.map(_.reproduce);
        }

        // return the new migrant pool
        working_groups.fold(empty_pop)(_ + _);
    }

    // run the simulation and display statistics
    def run: Unit = {
        itterate
        draw_graphs
    }

    // draw graphs and save as right.png and left.png
    def draw_graphs;

    // return a new empty population
    protected def empty_pop: PopType;

    // returns initial migrant pool
    protected def initialise: PopType;

    // assign individuals from the migrant population to groups
    // returns each group as a population in an Iterable
    protected def assign_to_groups(migrant_pool: PopType): Iterable[PopType];
}

