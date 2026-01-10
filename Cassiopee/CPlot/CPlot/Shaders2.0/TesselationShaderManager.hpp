/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _CPLOT_TESSELATIONSHADERMANAGER_HPP_
#define _CPLOT_TESSELATIONSHADERMANAGER_HPP_
#include <map>
#include <vector>
#include <utility>
#include <memory>
#include <string>
#include <iostream>
#include "TesselationControlShader.hpp"
#include "TesselationEvaluationShader.hpp"

namespace CPlot
{
  class TesselationShaderManager
  {
  public:
    TesselationShaderManager();
    ~TesselationShaderManager();
    
    /* load les shaders */
    int load();

    /**
     * @brief Add a tesselation shader for ho elements
     * @details Add a tesselation shader. A tesselation shader is composed with two shaders, one for control, one
     *          for evaluation.
     * 
     * @param tesselationBaseNameFile The base name for control and evaluation tesselation files. By example,
     *        triangle_p2 will load triangle_p2.tcs and triangle_p2.tes for control and evaluation shaders.
     * @return bool to tell if loading succeed
     */
    bool addFromFile( const char* tesselationBaseNameFile );

    unsigned short numberOfShaders() const
    { return (unsigned short)_shaderList.size(); }

    const std::pair<std::shared_ptr<TesselationControlShader>,
                    std::shared_ptr<TesselationEvaluationShader>>& activate(unsigned short id)
    {
      _currentActiveShader = id;
      return (*this)[id];
    }

    void deactivate()
    {
      _currentActiveShader = 0;
    }

    unsigned short currentShader() const
    { return _currentActiveShader; }
    
    bool is_activate() const { return _currentActiveShader > 0; }

    const std::pair<std::shared_ptr<TesselationControlShader>,
                    std::shared_ptr<TesselationEvaluationShader>>& operator [](unsigned short id) const 
    { return _shaderList[id-1]; }
  private:
    TesselationShaderManager(const TesselationShaderManager& shadMan);
    TesselationShaderManager& operator = (const TesselationShaderManager& shadMan);
    
    std::vector<std::pair<std::shared_ptr<TesselationControlShader>,std::shared_ptr<TesselationEvaluationShader>>> _shaderList;
    unsigned short _currentActiveShader;
  };
}
#endif
