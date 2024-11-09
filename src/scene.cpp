/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    v1 - Dec 01 2020
    Copyright (c) 2020 by Adrian Jarabo

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/scene.h>
#include <nori/bitmap.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/camera.h>
#include <nori/emitter.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

Scene::Scene(const PropertyList &) {
    m_accel = new Accel();
    m_enviromentalEmitter = 0;
}

Scene::~Scene() {

    delete m_accel;
    delete m_sampler;
    delete m_camera;
    delete m_integrator;

    for (size_t i=0; i<sss_accelerators.size(); ++i)
        delete sss_accelerators[i];

    for (size_t i=0; i<m_meshes.size(); ++i)
        delete m_meshes[i];
}

void Scene::activate() 
{
    // Build the discrete distribution based on radiance
    // m_emitterPDF.clear();
    // m_emitterPDF.reserve(m_meshes.size()+1);

    // Check if there's emitters attached to meshes, and
    // add them to the scene. 
    for(unsigned int i=0; i<m_meshes.size(); ++i )
    {
        if (m_meshes[i]->isEmitter())
        {
            m_meshes[i]->getEmitter()->setMesh(m_meshes[i]);
            m_emitters.push_back(m_meshes[i]->getEmitter());
            m_emitterPDF.append(m_meshes[i]->getEmitter()->getRadiance().maxCoeff());
        }

        if (m_meshes[i]->hasSubsurfaceScattering())
        {
            sss_meshes.push_back(m_meshes[i]);
        }
    }

    m_emitterPDF.normalize();

    m_accel->build();

    for (size_t i=0; i<sss_meshes.size(); ++i)
    {
        sss_accelerators.push_back(new Accel());
        sss_accelerators.back()->addMesh(sss_meshes[i]);
        sss_accelerators.back()->build();
    }

    if (!m_integrator)
        throw NoriException("No integrator was specified!");
    if (!m_camera)
        throw NoriException("No camera was specified!");
    
    if (!m_sampler) {
        /* Create a default (independent) sampler */
        m_sampler = static_cast<Sampler*>(
            NoriObjectFactory::createInstance("independent", PropertyList()));
    }

    if (m_emitters.size() == 0)
    {
        throw NoriException("No emitters were specified!");
    }

    cout << endl;
    cout << "Configuration: " << toString() << endl;
    cout << endl;
}


/// Sample emitter
// Emitter *Scene::sampleEmitter(Sampler* sampler, float &pdf) const
// {
//     float rnd = sampler->next1D();
//     float pdf2;
//     size_t index2 = m_emitterPDF.sample(rnd, pdf2);
//     cerr << "PDF m_emitterPDF: " << pdf2 << endl;
//     cerr << "m_emitterPDF.size: " << m_emitterPDF.size() << endl;
//     cerr << "m_emitterPDF.getSum: " << m_emitterPDF.getSum() << endl;
//     cerr << "m_emitterPDF.getNormalization: " << m_emitterPDF.getNormalization() << endl;


// 	auto const & n = m_emitters.size();
// 	size_t index = std::min(static_cast<size_t>(std::floor(n*rnd)), n - 1);
// 	pdf = 1. / float(n);
//     cerr << "PDF: " << pdf << endl;
// 	return m_emitters[index];
// }

void Scene::preprocess() {
    for (auto mesh : m_meshes)
    {
        mesh->getBSDF()->preprocess(m_sampler);
    }
}

Emitter *Scene::sampleEmitter(Sampler* sampler, float &pdf) const
{
    float rnd = sampler->next1D(); // Get a random value between 0 and 1

    // Sample from the PDF
    //size_t index = m_emitterPDF.sample(rnd, pdf); // Multi importance lighting
    size_t index = std::min(static_cast<size_t>(std::floor(m_emitters.size()*rnd)), m_emitters.size() - 1); // random light
    pdf = 1. / float(m_emitters.size());

    // Return the selected emitter
    return m_emitters[index];
}

float Scene::pdfEmitter(const Emitter *em) const {
    return 1. / float(m_emitters.size());
}


void Scene::addChild(NoriObject *obj, const std::string& name) {
    switch (obj->getClassType()) {
        case EMesh: {
                Mesh *mesh = static_cast<Mesh *>(obj);
                m_accel->addMesh(mesh);
                m_meshes.push_back(mesh);
            }
            break;
        
        case EEmitter: {
				Emitter *emitter = static_cast<Emitter *>(obj);
				if (emitter->getEmitterType() == EmitterType::EMITTER_ENVIRONMENT)
				{
					if (m_enviromentalEmitter)
						throw NoriException("There can only be one enviromental emitter per scene!");
					m_enviromentalEmitter = emitter;
				}
				
                m_emitters.push_back(emitter);
                m_emitterPDF.append(emitter->getRadiance().maxCoeff());
			}
            break;

        case ESampler:
            if (m_sampler)
                throw NoriException("There can only be one sampler per scene!");
            m_sampler = static_cast<Sampler *>(obj);
            break;

        case ECamera:
            if (m_camera)
                throw NoriException("There can only be one camera per scene!");
            m_camera = static_cast<Camera *>(obj);
            break;
        
        case EIntegrator:
            if (m_integrator)
                throw NoriException("There can only be one integrator per scene!");
            m_integrator = static_cast<Integrator *>(obj);
            break;

        default:
            throw NoriException("Scene::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
    }
}

Color3f Scene::getBackground(const Ray3f& ray) const
{
    if (!m_enviromentalEmitter)
        return Color3f(0);

    EmitterQueryRecord lRec(m_enviromentalEmitter, ray.o, ray.o + ray.d, Normal3f(0, 0, 1), Vector2f());
	return m_enviromentalEmitter->eval(lRec);
}


std::string Scene::toString() const 
{
    std::string meshes;
    for (size_t i=0; i<m_meshes.size(); ++i) {
        meshes += std::string("  ") + indent(m_meshes[i]->toString(), 2);
        if (i + 1 < m_meshes.size())
            meshes += ",";
        meshes += "\n";
    }

	std::string lights;
	for (size_t i = 0; i < m_emitters.size(); ++i) {
		lights += std::string("  ") + indent(m_emitters[i]->toString(), 2);
		if (i + 1 < m_emitters.size())
			lights += ",";
		lights += "\n";
	}


    return tfm::format(
        "Scene[\n"
        "  integrator = %s,\n"
        "  sampler = %s\n"
        "  camera = %s,\n"
        "  meshes = {\n"
        "  %s  }\n"
		"  emitters = {\n"
		"  %s  }\n"
        "]",
        indent(m_integrator->toString()),
        indent(m_sampler->toString()),
        indent(m_camera->toString()),
        indent(meshes, 2),
		indent(lights, 2)
    );
}

NORI_REGISTER_CLASS(Scene, "scene");
NORI_NAMESPACE_END
