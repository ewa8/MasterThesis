
/*
 *  nestml_7572c325c24847f4afc90a8b52d50e81_module.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  2022-11-16 08:28:56.405724
 */

#ifndef NESTML_7572C325C24847F4AFC90A8B52D50E81_MODULE_H
#define NESTML_7572C325C24847F4AFC90A8B52D50E81_MODULE_H

#include "slimodule.h"
#include "slifunction.h"

#include "nest.h"
#include "nest_impl.h"


/**
* Class defining your model.
* @note For each model, you must define one such class, with a unique name.
*/
class nestml_7572c325c24847f4afc90a8b52d50e81_module : public SLIModule
{
public:
  // Interface functions ------------------------------------------

  /**
   * @note The constructor registers the module with the dynamic loader.
   *       Initialization proper is performed by the init() method.
   */
  nestml_7572c325c24847f4afc90a8b52d50e81_module();

  /**
   * @note The destructor does not do much in modules.
   */
  ~nestml_7572c325c24847f4afc90a8b52d50e81_module();

  /**
   * Initialize module by registering models with the network.
   * @param SLIInterpreter* SLI interpreter
   */
  void init( SLIInterpreter* );

  /**
   * Return the name of your model.
   */
  const std::string name() const;

public:
  // Classes implementing your functions -----------------------------

};

#endif