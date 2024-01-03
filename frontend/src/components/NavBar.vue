<template>
  <v-app-bar :color="$route.meta.navBarColor" class="nav-bar" :elevation="4" :style="{ transition: 'background-color 0.35s' }">
    <v-app-bar-nav-icon icon="mdi-react" @click="redirectHome"></v-app-bar-nav-icon>
    <v-app-bar-title class="nav-bar-title" @click="redirectHome">Reaction Playground</v-app-bar-title>

    <!--    Add the nav bar route buttons-->
    <v-btn v-for="item in navBarItems" :key="item.text" :prepend-icon="item.icon" variant="text" :to="{ name: item.routeName }">
      {{ item.text }}
    </v-btn>
    <v-divider
        class="mx-3 align-self-center"
        length="24"
        thickness="2"
        vertical
    ></v-divider>
    <v-btn icon>
      <v-icon>mdi-cog</v-icon>
    </v-btn>
    <v-btn icon>
      <v-icon>mdi-help-circle</v-icon>
    </v-btn>

  </v-app-bar>
</template>

<script>
import { defineComponent } from 'vue'

export default defineComponent({
  name: 'NavBar',
  data () {
    return {
      navBarItems: [
        { icon: 'mdi-flask-outline', text: 'Playground', routeName: 'playgroundMode' },
        { icon: 'mdi-puzzle-outline', text: 'Solver', routeName: 'solverMode' },
        { icon: 'mdi-text-box-search-outline', text: 'Classifier', routeName: 'molClassifier' }
      ]
    }
  },
  methods: {
    redirectHome () {
      this.$router.push({ name: 'home' })
    }
  }
})
</script>

<style>
.v-app-bar-title {
  cursor: pointer;
}
.v-btn--active {
  background-color: transparent;
  color: inherit;
  border-bottom: 2px solid white; /* Set the color and style of the underline */
}
</style>
