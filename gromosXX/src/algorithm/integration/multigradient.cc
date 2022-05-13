/**
 * @file multigradient.cc
 * multigradient implementation
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"


#include "multigradient.h"



#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

algorithm::Multi_Gradient::ControlPoint operator+(
        const algorithm::Multi_Gradient::ControlPoint & lhs,
        const algorithm::Multi_Gradient::ControlPoint & rhs) {
  return algorithm::Multi_Gradient::ControlPoint(lhs.time + rhs.time, lhs.value + rhs.value);
}

algorithm::Multi_Gradient::ControlPoint operator*(
        const algorithm::Multi_Gradient::ControlPoint & lhs, double rhs) {
  return algorithm::Multi_Gradient::ControlPoint(lhs.time * rhs, lhs.value * rhs);
}

void algorithm::Multi_Gradient::LinearInterpolation::echo(std::ostream & os) const {
  os << "    Linear interpolation using " << control_points().size() << " control points: " << std::endl;
  for (unsigned int j = 0; j < control_points().size(); ++j) {
    os.precision(4);
    os << "    - " << std::setw(8) << control_points()[j].time << ": "
            << std::setw(8) << control_points()[j].value << std::endl;
  }
}

double algorithm::Multi_Gradient::LinearInterpolation::get_value(double time) const {

  if (time < control_points()[0].time)
    return control_points()[0].value;

  if (time >= control_points()[control_points().size()-1].time)
    return control_points()[control_points().size()-1].value;

  for(unsigned int i = 0; i < control_points().size()-1; ++i) {
    if (control_points()[i].time <= time && control_points()[i + 1].time > time) {
      double t = time - control_points()[i].time;
      double dt =  control_points()[i + 1].time - control_points()[i].time;
      double dv = control_points()[i + 1].value - control_points()[i].value;
      return control_points()[i].value + t*dv/dt;
    }
  }
  return control_points()[control_points().size()-1].value;;
}

void algorithm::Multi_Gradient::SplineInterpolation::init() {
  // get arrays for x/y values
  for (unsigned int j = 0; j < control_points().size(); ++j) {
    times.push_back(control_points()[j].time);
    values.push_back(control_points()[j].value);
  }
  // allocate and initalize the splines
  spline = gsl_spline_alloc(gsl_interp_cspline, control_points().size());
  gsl_spline_init(spline, &times[0], &values[0], control_points().size());
}

void algorithm::Multi_Gradient::SplineInterpolation::echo(std::ostream & os) const {
  os << "    Cubic spline interpolation using " << control_points().size() << " control points: " << std::endl;
  for (unsigned int j = 0; j < control_points().size(); ++j) {
    os.precision(4);
    os << "    - " << std::setw(8) << control_points()[j].time << ": "
            << std::setw(8) << control_points()[j].value << std::endl;
  }
}

double algorithm::Multi_Gradient::SplineInterpolation::get_value(double time) const {

  if (time < control_points()[0].time)
    return control_points()[0].value;

  if (time >= control_points()[control_points().size()-1].time)
    return control_points()[control_points().size()-1].value;

  return gsl_spline_eval(spline, time, acc);
}

double algorithm::Multi_Gradient::Bezier::get_value(double time) const {

  if (time < control_points()[0].time)
    return control_points()[0].value;

  if (time > control_points()[control_points().size()-1].time)
    return control_points()[control_points().size()-1].value;

  // de Casteljau algorithm
  double t = (time - control_points()[0].time) / (control_points()[control_points().size()-1].time - control_points()[0].time);
  double one_minus_t = 1.0 - t;

  std::vector<algorithm::Multi_Gradient::ControlPoint> q = control_points();
  for(unsigned int k = 1; k < q.size(); ++k) {
    for(unsigned int i = 0; i < q.size() - k; ++i) {
      q[i] = q[i]*one_minus_t + q[i+1] * t;
    }
  }
  return q[0].value;
}

void algorithm::Multi_Gradient::Bezier::echo(std::ostream & os) const {
  os << "    Bezier using " << control_points().size() << " control points: " << std::endl;
  for (unsigned int j = 0; j < control_points().size(); ++j) {
    os.precision(4);
    os << "    - " << std::setw(8) << control_points()[j].time << ": "
            << std::setw(8) << control_points()[j].value << std::endl;
  }
}

double algorithm::Multi_Gradient::Oscillation::get_value(double time) const {
  const double & A =  control_points()[0].time;
  const double & T = control_points()[0].value;
  const double & DeltaT = control_points()[1].time;
  const double & b = control_points()[1].value;
  return A * sin(2.0 * math::Pi * (time - DeltaT) / T) + b;
}

void algorithm::Multi_Gradient::Oscillation::echo(std::ostream & os) const {
  os << "    Oscillation A sin[2Pi/T (t - dt)] + b" << std::endl;
  os << "    - A:  " << std::setw(8) << control_points()[0].time << std::endl;
  os << "    - T:  " << std::setw(8) << control_points()[0].value << std::endl;
  os << "    - dt:  " << std::setw(8) << control_points()[1].time << std::endl;
  os << "    - b:  " << std::setw(8) << control_points()[1].value << std::endl;
}

std::string algorithm::Multi_Gradient::Curve::plot_ascii(double start_time, double end_time,
        unsigned int width, unsigned int height,
        std::string x_label, std::string y_label, const std::string & indent) const {
  char matrix[width][height];
  for(unsigned int y = 0; y < height; ++y) {
    for(unsigned int x = 0; x < width; ++x) {
      matrix[x][y] = ' ';
    }
  }

  double dt = (end_time - start_time) / width;
  double values[width];
  for(unsigned int i = 0; i < width; ++i) {
    double time = start_time + dt * i;
    values[i] = get_value(time);
  }

  double min_val = values[0];
  double max_val = values[0];
  for(unsigned int i = 0; i < width; ++i) {
    min_val = std::min(min_val, values[i]);
    max_val = std::max(max_val, values[i]);
  }

  double spread = max_val - min_val;

  for(unsigned int x = 0; x < width; ++x) {
    int y = 0;
    if (spread)
      y = int(height * (values[x] - min_val) / spread);
    else
      y = 0;

    y = std::max(0, y);
    y = std::min(int(height-1), y);

    matrix[x][y] = '*';
  }

  std::ostringstream os;
  unsigned int i = 0;
  y_label = y_label.substr(0, height);
  x_label = x_label.substr(0, width);
  unsigned int lab_offset = (height - y_label.size())/2;
  for(int y = height - 1; y >= 0; --y,++i) {
    os << indent;

    if (i >= lab_offset && (i-lab_offset) < y_label.size())
      os << y_label[i-lab_offset];
    else
      os << " ";

    os << " ";

    if (i == 0 && spread) {
      os.precision(4);
      os.setf(std::ios::right);
      os << std::setw(10) << std::scientific << max_val << "|";
    } else if (y == 0) {
      os.precision(4);
      os.setf(std::ios::right);
      os << std::setw(10) << std::scientific << min_val << "|";
    } else {
      os << "          | ";
    }

    for(unsigned int x = 0; x < width; ++x) {
      os << matrix[x][y];
    }
    os << std::endl;
  }
  os << indent << "            +-";
  for(unsigned int x = 0; x < width-1; ++x)
    os << "-";
  os << ">" << std::endl;

  os << indent << "            ";
  os.precision(4);
  os.setf(std::ios::left);
  os << std::setw(8) << start_time;
  for(unsigned int i = 0; i < width - 16; ++i) {
    os << ' ';
  }
  os.setf(std::ios::right);
  os << std::setw(8) << end_time << std::endl;

  if (x_label.size()) {
    os << indent << "          ";
    for(unsigned int i = 0; i < (width - x_label.size()) / 2; ++i) {
      os << ' ';
    }
    os << x_label << std::endl;
  }
  return os.str();
}

bool parse_var(const std::string & str, std::string & var, int & index) {
  std::string::size_type bra = str.find_first_of("["),
          ket = str.find_first_of("]");
  if (bra == std::string::npos && ket == std::string::npos) {
    var = str; index = -1;
    return true;
  }

  if (ket != str.size() - 1)
    return false;

  var = str.substr(0, bra);
  std::string np = str.substr(bra+1, str.size() - bra - 2);
  std::istringstream num_part(np);
  if (num_part >> index) {
    return true;
  }

  return false;

}
bool check_variable(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation &sim,
        const std::string & var) {
  std::string name;
  int index = -1;
  if (!parse_var(var, name, index)) {
    io::messages.add("Malformed variable " + var,
                "Multi_Gradient", io::message::error);
    return false;
  }

  if (name == "TEMP0") {
    --index;
    if (index < 0 || index > int(sim.multibath().size())) {
      std::ostringstream msg;
      msg << "There is no bath with index " << (index + 1);
      io::messages.add(msg.str(), "Multi_Gradient", io::message::error);
      return false;
    }
    if (!sim.param().multibath.couple) {
      io::messages.add("Changing for TEMP0 without temperature coupling makes no sense",
              "Multi_Gradient", io::message::error);
      return false;
    }
    return true;
  } else if (name == "CPOR") {
    if (sim.param().posrest.posrest == simulation::posrest_off ||
        sim.param().posrest.posrest == simulation::posrest_const) {
      io::messages.add("Changing for CPOR without position restraining makes no sense",
              "Multi_Gradient", io::message::error);
      return false;
    }
    return true;
  } else if (name == "CDIR") {
    if (!sim.param().distanceres.distanceres) {
      io::messages.add("Changing for CDIR without position restraining makes no sense",
              "Multi_Gradient", io::message::error);
      return false;
    }
    return true;
  } else if (name == "RESO") {
    if (sim.param().xrayrest.xrayrest == simulation::xrayrest_off) {
      io::messages.add("Changing for RESO without position restraining makes no sense",
              "Multi_Gradient", io::message::error);
      return false;
    }
    return true;
  } else if (name == "CXR") {
    if (sim.param().xrayrest.xrayrest == simulation::xrayrest_off) {
      io::messages.add("Changing for CXR without X-ray restraining makes no sense",
              "Multi_Gradient", io::message::error);
      return false;
    }
    return true;
  } else if (name == "COPR") {
    if (sim.param().orderparamrest.orderparamrest == simulation::oparam_restr_off) {
      io::messages.add("Changing for COPR without order-parameter restraining makes no sense",
              "Multi_Gradient", io::message::error);
      return false;
    }
    return true;
  }

  io::messages.add("Unkown variable " + var,
                "Multi_Gradient", io::message::error);
  return false;
}

int algorithm::Multi_Gradient
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation &sim)
{

  m_timer.start();
  // loop over curves
  for(std::vector<algorithm::Multi_Gradient::Curve*>::const_iterator it = curves.begin(),
            to = curves.end(); it != to; ++it) {
    std::string name;
    int index = 0;
    parse_var((*it)->variable, name, index);

    if (name == "TEMP0") {
      sim.multibath().bath(index-1).temperature = (*it)->get_value(sim.time());
    } else if (name == "CPOR") {
      sim.param().posrest.force_constant = (*it)->get_value(sim.time());
    } else if (name == "CDIR") {
      sim.param().distanceres.K = (*it)->get_value(sim.time());
    } else if (name == "RESO") {
      sim.param().xrayrest.resolution = (*it)->get_value(sim.time());
    } else if (name == "CXR") {
      sim.param().xrayrest.force_constant = (*it)->get_value(sim.time());
    } else if (name == "COPR") {
      sim.param().orderparamrest.K = (*it)->get_value(sim.time());
    }
  }
  m_timer.stop();

  return 0;
}

int algorithm::Multi_Gradient
::init(topology::Topology &topo,
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os,
		     bool quiet)
{

  for(unsigned int i = 0; i < sim.param().multigradient.variable.size(); ++i) {
    algorithm::Multi_Gradient::Curve * c = nullptr;
    
    // check
    if (!check_variable(topo, conf, sim, sim.param().multigradient.variable[i])) {
        return 1;
    }

    switch(sim.param().multigradient.functional_form[i]) {
      case 0:
        c = new algorithm::Multi_Gradient::LinearInterpolation;
        break;
      case 1:
        c = new algorithm::Multi_Gradient::SplineInterpolation;
        break;
      case 2:
        c = new algorithm::Multi_Gradient::Bezier;
        break;
      case 3:
        c = new algorithm::Multi_Gradient::Oscillation;
        break;
      default:
        io::messages.add("Unkown functional form",
                "Multi_Gradient", io::message::error);
        return 1;
    }

    c->variable = sim.param().multigradient.variable[i];

    for(unsigned int j = 0; j < sim.param().multigradient.control_points[i].size(); ++j) {
      algorithm::Multi_Gradient::ControlPoint p(sim.param().multigradient.control_points[i][j].first,
              sim.param().multigradient.control_points[i][j].second);
      c->add_control_point(p);
    }

    // do not check for oscillation
    if (sim.param().multigradient.functional_form[i] <= 2) {
      double time = c->control_points()[0].time;
      for (unsigned int j = 1; j < c->control_points().size(); ++j) {
        if (c->control_points()[j].time < time) {
          io::messages.add("Control points have to be successive in time.",
                  "Multi_Gradient", io::message::error);
          return 1;
        }
        time = c->control_points()[j].time;
      }
    } else {
      if (c->control_points().size() != 2) {
        io::messages.add("An osicillation takes 4 parameters (i.e. 2 control points!)",
                  "Multi_Gradient", io::message::error);
          return 1;
      }
    }
    // init the curve
    c->init();
    curves.push_back(c);
  }

  if (!quiet) {
    os << "MULTIGRADIENT" << std::endl;
    os << "Multigradient is " <<
            (sim.param().multigradient.multigradient ? "enabled" : "disabled") << "." << std::endl;
    os << "Number of curves: " << curves.size() << std::endl;
    unsigned int i = 1;
    for(std::vector<algorithm::Multi_Gradient::Curve*>::const_iterator it = curves.begin(),
            to = curves.end(); it != to; ++it, ++i) {
      os << "Curve " << i << " affecting variable " << (*it)->variable << std::endl;
      (*it)->echo(os);

      double start_time = sim.time();
      double end_time = start_time + sim.param().step.number_of_steps * sim.time_step_size();
      if (sim.param().multigradient.print_graph) {
        os << std::endl;
        os << (*it)->plot_ascii(start_time, end_time, 50, 10, "TIME", (*it)->variable, "        ");
      }

      if (sim.param().multigradient.print_curve) {
        os << "    100 exact values: " << std::endl;
        os.precision(8);
        os << std::scientific;
        for(unsigned int i = 0; i < 100; ++i) {
          double time = start_time + i * (end_time - start_time) / 100;
          os << "      " << std::setw(15) << time << std::setw(15) << (*it)->get_value(time) << std::endl;
        }
      }

      os << std::endl;
    }
    os << "END" << std::endl;
  }
  return 0;
}

algorithm::Multi_Gradient::~Multi_Gradient() {
  for(std::vector<algorithm::Multi_Gradient::Curve*>::const_iterator it = curves.begin(),
            to = curves.end(); it != to; ++it)
    delete *it;
}
