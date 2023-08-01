// Compile using:
//
// clang++ -O3 -Werror -Wall -std=c++20 lunar.cpp

/*
Can we compute,
    "if I did a suicide burn, how much fuel would that take?" If that's more
than we have, then we can end here.;
*/

#include <cassert>
#include <cmath>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace std;

constexpr double EMPTY_WEIGHT = 16500.0;         // N, lbs
constexpr double INITIAL_TOTAL_WEIGHT = 32500.0; // M, lbs
constexpr double INITIAL_ALTITUDE = 120;         // A, miles
constexpr double INITIAL_VELOCITY = 1.0;         // V, miles/sec

constexpr double GRAVITY = 0.001;        // G, miles / second^2
constexpr double EXHAUST_VELOCITY = 1.8; // Z, miles / second.

constexpr double INVALID_VELOCITY = 1e6;

void print_vector(const vector<double> &vec) {
  printf("{");
  for (const double &elem : vec) {
    printf("%.0f", elem);
    if (&elem != &vec.back()) {
      printf(", ");
    }
  }
  printf("}");
}

struct Lander {
  double time;         // Total elapsed time of simulation (sec)
  double total_weight; // Dry weight plus fuel (lbs)
  double empty_weight; // Dry weight (lbs)
  double altitude;     // Miles above surface.
  double
      velocity; // Downward velocity, so opposite sign of altitude (miles/sec)
  double velocity_when_out_of_fuel{INVALID_VELOCITY};
  bool on_moon{false};

  Lander()
      : time(0), total_weight(INITIAL_TOTAL_WEIGHT), empty_weight(EMPTY_WEIGHT),
        altitude(INITIAL_ALTITUDE), velocity(INITIAL_VELOCITY) {}

  double fuel_weight() const { return total_weight - empty_weight; }

  bool fuel_empty() const { return fuel_weight() < 0.001; }

  double velocity_mph() const { return 3600 * velocity; }

  bool ran_out_of_fuel() const {
    return velocity_when_out_of_fuel != INVALID_VELOCITY;
  }

  bool out_and_up() const {
    return ran_out_of_fuel() && velocity_when_out_of_fuel < 0;
  }

  void set_on_moon() {
    on_moon = true;
    altitude = 0;
  }

  void update(double elapsed, double new_altitude, double new_velocity,
              double fuel_rate) {
    time += elapsed;
    total_weight -= elapsed * fuel_rate;
    altitude = new_altitude;
    velocity = new_velocity;
  }

  void update_when_out_of_fuel() {
    velocity_when_out_of_fuel = velocity;
    double elapsed =
        (sqrt(velocity * velocity + 2 * altitude * GRAVITY) - velocity) /
        GRAVITY;
    velocity += GRAVITY * elapsed;
    time += elapsed;
    set_on_moon();
  }

  // Note: doesn't detect crash in the middle of timestep.  Caller needs to do
  // that.
  pair<double, double> apply_thrust(double elapsed, double fuel_rate) const {
    // Fraction of total weight reduced.
    double q = elapsed * fuel_rate / total_weight;
    double q2 = q * q;
    double q3 = q2 * q;
    double q4 = q3 * q;
    double q5 = q4 * q;

    // Return these two values.
    double end_velocity = velocity + GRAVITY * elapsed +
                          // Tsiolkovsky rocket equation
                          EXHAUST_VELOCITY *
                              // Taylor series for ln(1 - q)
                              (-q - q2 / 2 - q3 / 3 - q4 / 4 - q5 / 5);

    double end_altitude =
        altitude - GRAVITY * elapsed * elapsed / 2 - velocity * elapsed +
        EXHAUST_VELOCITY * elapsed *
            // Integral of above, i.e. integral of the rocket equation.
            (q / 2 + q2 / 6 + q3 / 12 + q4 / 20 + q5 / 30);

    return {end_altitude, end_velocity};
  }

  void loop_until_on_the_moon(double elapsed, double fuel_rate,
                              double &turn_time) {
    while (elapsed >= 0.005) {
      elapsed =
          2 * altitude /
          (velocity +
           sqrt(velocity * velocity +
                2 * altitude *
                    (GRAVITY - EXHAUST_VELOCITY * fuel_rate / total_weight)));
      double new_altitude, new_velocity;
      tie(new_altitude, new_velocity) = apply_thrust(elapsed, fuel_rate);
      turn_time -= elapsed;
      update(elapsed, new_altitude, new_velocity, fuel_rate);
    }
    set_on_moon();
  }

  void print_header() const {
    puts("TIME,SECS   ALTITUDE,MILES+FEET   VELOCITY,MPH   FUEL,LBS   FUEL "
         "RATE");
  }

  void print(double fuel_rate) const {
    printf("%7.0f%16.0f%7.0f%15.2f%12.1f%10.1f\n", time, trunc(altitude),
           5280 * (altitude - trunc(altitude)), velocity_mph(), fuel_weight(),
           fuel_rate);
  }
  void print_on_moon() const {
    if (ran_out_of_fuel()) {
      printf("Ran out of fuel with velocity %f\n", velocity_when_out_of_fuel);
    }
    printf("On the moon at %f sec\n", time);
    printf("Impact velocity of %f m.p.h.\n", velocity_mph());
    printf("Fuel left %f lbs\n", fuel_weight());
  }
};

// "true" means on the moon.
bool do_turn(const double fuel_rate, Lander &lander) {
  double turn_time = 10.0;
  for (;;) {
    if (lander.fuel_empty()) {
      // printf("Out of fuel.\n");
      lander.update_when_out_of_fuel();
      return true;
    }

    if (turn_time < 0.001) {
      return false;
    }

    double elapsed = turn_time;

    // Are we going to run out of fuel before elapsed?
    if (elapsed * fuel_rate > lander.fuel_weight()) {
      elapsed = lander.fuel_weight() / fuel_rate;
    }

    double new_altitude, new_velocity;
    tie(new_altitude, new_velocity) = lander.apply_thrust(elapsed, fuel_rate);

    if (new_altitude <= 0) {
      lander.loop_until_on_the_moon(elapsed, fuel_rate, turn_time);
      return true;
    }

    if (lander.velocity > 0 && new_velocity < 0) {
      // We were heading down, are now heading up.  Did we collide with the
      // lunar surface?
      while (true) {
        // I'm not sure where this formula comes from.  On the one hand, it
        // shows a lot of knowledge of the physics and math of the
        // situation, e.g. the units all work out and it's a straight
        // forward expression.  On the other hand, the simplest way to
        // derive the "elapsed" where our altitude is lowest and velocity is
        // zero leads to a significantly more accurate expression.

        // This should be a method of Lander.
        const double ratio = lander.total_weight / fuel_rate; // sec
        const double W = (1 - ratio * GRAVITY / EXHAUST_VELOCITY) / 2;
        const double velocity_prime = lander.velocity / EXHAUST_VELOCITY;
        elapsed =
            ratio * velocity_prime / (W + sqrt(W * W + velocity_prime)) + 0.05;

        tie(new_altitude, new_velocity) =
            lander.apply_thrust(elapsed, fuel_rate);

        if (new_altitude <= 0) {
          lander.loop_until_on_the_moon(elapsed, fuel_rate, turn_time);
          return true;
        }
        turn_time -= elapsed;
        lander.update(elapsed, new_altitude, new_velocity, fuel_rate);
        if (new_velocity > 0 || lander.velocity <= 0) {
          break;
        }
      }
    } else {
      turn_time -= elapsed;
      lander.update(elapsed, new_altitude, new_velocity, fuel_rate);
    }
  }
}

struct Result {
  bool on_moon;
  double velocity_mph;
  double fuel;
  double time;
  bool ran_out_of_fuel;
};

Lander do_run_up_to(const vector<double> &fuel_schedule, size_t num_iter) {
  Lander lander;

  for (size_t i = 0; i < num_iter; ++i) {
    const double fuel_rate = fuel_schedule[i];
    bool on_moon = do_turn(fuel_rate, lander);
    if (on_moon) {
      return lander;
    }
  }
  return lander;
}

// TODO: Return Lander, and we can get rid of the Result struct as well as
// do_run_up_to.
Result do_run(const vector<double> &fuel_schedule, ssize_t max_iter = -1,
              bool verbose = false) {
  Lander lander;

  if (verbose) {
    lander.print_header();
  }

  if (max_iter < 0 || max_iter > (ssize_t)fuel_schedule.size()) {
    max_iter = fuel_schedule.size();
  }

  for (ssize_t i = 0; i < max_iter; ++i) {
    const double fuel_rate = fuel_schedule[i];
    if (verbose) {
      lander.print(fuel_rate);
    }
    bool on_moon = do_turn(fuel_rate, lander);
    if (on_moon) {
      // printf("ON MOON.\n");
      if (verbose) {
        lander.print_on_moon();
      }
      return Result{true, lander.velocity_mph(), lander.fuel_weight(),
                    lander.time, lander.ran_out_of_fuel()};
    } else {
      // printf("not on moon.\n");
    }
  }

  // What to do when schedule runs out?
  return Result{false, lander.velocity_mph(), lander.fuel_weight(), lander.time,
                lander.ran_out_of_fuel()};
}

double last_fuel_brute_force(vector<double> &fuel_schedule) {
  double max_fuel = -1;
  double max_fuel_rate = -1;

  const Lander lander{do_run_up_to(fuel_schedule, fuel_schedule.size() - 1)};

  if (lander.on_moon) {
    fuel_schedule.back() = 0;
    if (lander.velocity_mph() > 1) {
      return lander.fuel_weight();
    } else {
      return -1;
    }
  }

  for (int fuel_rate = 0; fuel_rate <= 200; fuel_rate++) {
    if (fuel_rate == 1) {
      fuel_rate = 8;
    }

    Lander new_lander{lander};
    do_turn(fuel_rate, new_lander);

    if (!new_lander.on_moon || new_lander.velocity_mph() > 1) {
      continue;
    }

    if (new_lander.fuel_weight() > max_fuel) {
      max_fuel = new_lander.fuel_weight();
      max_fuel_rate = fuel_rate;
    }
  }

  fuel_schedule.back() = max_fuel_rate;

  return max_fuel;
}

double nth_last_brute_force(vector<double> &fuel_schedule, int n) {
  assert(n > 1);

  double max_fuel = -1;
  vector<double> max_fuel_schedule;

  for (int fuel_rate = 0; fuel_rate <= 200; fuel_rate++) {
    if (fuel_rate == 1) {
      fuel_rate = 8;
    }

    fuel_schedule[fuel_schedule.size() - n] = fuel_rate;

    for (int i = fuel_schedule.size() - (n - 1);
         i < (ssize_t)fuel_schedule.size(); i++) {
      fuel_schedule[i] = 200;
    }

    Result result = do_run(fuel_schedule);
    if (!result.ran_out_of_fuel && result.on_moon && result.velocity_mph > 1) {
      continue;
    }

    if (n >= 6) {
      for (int i = 0; i < 10 - n; ++i) {
        printf("  ");
      }
      printf("n: %d, trying %d\n", n, fuel_rate);
    }

    double fuel;
    if (n == 2) {
      fuel = last_fuel_brute_force(fuel_schedule);
    } else {
      fuel = nth_last_brute_force(fuel_schedule, n - 1);
    }

    if (fuel > max_fuel) {
      if (n >= 7) {
        printf("***** n=%d New max fuel: %f @ fuel rate %d\n", n, fuel,
               fuel_rate);
        print_vector(fuel_schedule);
        printf("\n");
      }
      max_fuel = fuel;
      max_fuel_schedule = fuel_schedule;
    }
  }

  if (max_fuel >= 0) {
    fuel_schedule = max_fuel_schedule;
  }

  return max_fuel;
}

vector<double> serge_schedule = {0,   0,   0,   0,   0,   0,   0,   200,
                                 200, 189, 189, 177, 183, 198, 198, 10};

// For 8th entry @ 70 sec:
//   164.3142678 too low.
//   164.3142679 works.
// vector<double> my_schedule = {0,   0,   0,   0,   0,   0,   0,   164.3142679,
//                              200, 200, 200, 200, 200, 200, 200, 0};

// Fuel left:  602.65447 lbs.
// vector<double> my_schedule = {0,   0,   0,   0,   0,   0,   0,   167,
//                              200, 200, 200, 200, 200, 200, 166, 11.6};

// Lands with 3.24 mph
// vector<double> my_schedule = {0,   0,   0,   0,   0,   0,   0,  166,
//                               200, 200, 200, 200, 200, 200, 176.72};

// I think the problem is, his solver for intersecting the surface in the middle
// of a turn is buggy, so you easily end up choosing between 3 mph and not
// landing at all.  So you can't land in 150 seconds, and have to go for that
// extra turn.

/*
vector<double> my_schedule = {0,   0,   0,   0,   0,   0,   0,   180,
                              200, 200, 200, 200, 200, 200, 166, 11};
*/

// Remaining fuel: 655.6
vector<double> better_schedule = {0,   0,   0,   0,   0,   0,   0,   200,
                                  200, 189, 189, 166, 199, 198, 192, 50};

// Remaining fuel: 658.458338
vector<double> best_schedule = {0,   0,   0,   0,   0,   0,   0,   190,
                                189, 191, 192, 199, 200, 195, 178, 62};

vector<double> my_schedule = {0,   0,   0, 0, 0, 0, 0, 200,
                              200, 200, 0, 0, 0, 0, 0};

// There is no perfect landing (i.e. final velocity < 1mph) in 15 steps.  At
// least, if the first 7 steps all use zero fuel.

// First of 176, 198 gives
// Remaining fuel: 658.436

// First of 177, 194 gives
// Remaining fuel: 658.441

// First of 179, 194 gives
// Remaining fuel: 658.447

// First of 180 gives 658.441

// First of 181 gives 658.445049

// First of 182, 199, 198 gives 658.452443

// First of 183, 198, 199 gives at least 658.457687

// First of 184 at least 658.452
// {0, 0, 0, 0, 0, 0, 0, 184, 196, 195, 199, 197, 179, 198, 186, 61}

// First of 185 gives 658.452

// First of 186 gives 658.455

// First of 187 gives 658.45196
// {0, 0, 0, 0, 0, 0, 0, 187, 200, 197, 187, 189, 195, 187, 192, 62}

// First of 188 gives 658.457484

// First of 189 gives 658.449360

// First of 190 gives 658.458338  *****
// {0, 0, 0, 0, 0, 0, 0, 190, 189, 191, 192, 199, 200, 195, 178, 62}

// First of 191 gives 658.453568

// First of 192 gives 658.453607

// First of 193 gives 658.457310

// First of 194 gives 658.457547

// First of 195 gives 658.458219  ***

// First of 196 gives 658.455872

// 197 gives 658.454999

// 198 gives 658.457991

// 199 gives 658.457615

// 200 gives 658.458152
// {0, 0, 0, 0, 0, 0, 0, 200, 199, 184, 185, 185, 194, 195, 192, 62}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("Usage: %s first-nonzero-fuel-rate\n", argv[0]);
    return 1;
  }

  double fuel_rate = stod(argv[1]);
  if (fuel_rate < 8 || fuel_rate > 200) {
    fprintf(stderr, "Fuel rate %f not allowed.\n", fuel_rate);
    return 1;
  }
  my_schedule[6] = fuel_rate;

  printf("*****  Default Fuel Schedule\n");
  do_run(my_schedule, -1, true);

  double fuel = nth_last_brute_force(my_schedule, 8);
  if (fuel >= 0) {
    printf("*****  Running best schedule, ending fuel should be %f\n", fuel);
    do_run(my_schedule, -1, true);
  } else {
    printf("*****  Couldn't find working schedule.\n");
  }

  return 0;
}
