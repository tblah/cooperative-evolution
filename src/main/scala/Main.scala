package CooperativeEvolution

object Main extends App {
  println("Cooperative Evolution")

  // based on example at https://github.com/scalanlp/breeze/wiki/Quickstart
  /*val figure = new Figure("Normal Distribution", 1, 1)
  val subplot = figure.subplot(0);
  val gaussian = breeze.stats.distributions.Gaussian(0, 1)
  subplot += hist(gaussian.sample(1000000), 100)*/

  val p = new Paper;
  p.run;
}
