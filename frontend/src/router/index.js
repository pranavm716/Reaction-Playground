import { createRouter, createWebHistory } from 'vue-router'
import HomeView from '../views/HomeView.vue'
import PlaygroundMode from '../views/PlaygroundMode.vue'
import SolverMode from '../views/SolverMode.vue'
import MolClassifier from '../views/MolClassifier.vue'

const routes = [
  {
    path: '/',
    name: 'home',
    component: HomeView,
    meta: { navBarColor: 'light-blue-darken-4' }
  },
  {
    path: '/playground-mode',
    name: 'playgroundMode',
    component: PlaygroundMode,
    meta: { navBarColor: 'light-blue-darken-3' }
  },
  {
    path: '/solver-mode',
    name: 'solverMode',
    component: SolverMode,
    meta: { navBarColor: 'red-darken-1' }
  },
  {
    path: '/analyze-molecule',
    name: 'molClassifier',
    component: MolClassifier,
    meta: { navBarColor: 'teal' }
  }
]

const router = createRouter({
  history: createWebHistory(process.env.BASE_URL),
  routes
})

export default router
